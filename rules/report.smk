def template(s, mapping):
  return re.sub("{{(.*?)}}", lambda m: str(mapping[m.group(1)]), s)


rule report_render:
  input: report=join_path("results/report/report.Rmd"),
         config=join_path("results/report/config.Rmd"),
         #read_summary=join_path("results/report/read_summary.Rmd"),
         feature_summary=join_path("results/report/feature_summary.Rmd"),
         session=join_path("results/report/session.Rmd"),
  output: join_path("results/report/report.html"),
  log: join_path("logs/report/render.log"),
  shell: """
    Rscript {workflow.basedir}/scripts/render_report.R --format html_document --output {output} {input.report}
  """


rule report_template:
  input: f"{workflow.basedir}/report/report.Rmd",
  output: join_path("results/report/report.Rmd")
  params:
    mapping={"OUTPUT": join_path("results/report/")},
  run:
    with open(input[0], "r") as in_file:
      s = in_file.read()
      d = template(s, params.mapping)
      with open(output[0], "w") as out_file:
        out_file.write(d)


rule report_config_template:
  input: f"{workflow.basedir}/report/config.Rmd",
  output: join_path("results/report/config.Rmd")
  params:
    mapping={"OUTPUT": join_path("results/report/")},
  run:
    with open(input[0], "r") as in_file:
      s = in_file.read()
      d = template(s, params.mapping)
      with open(output[0], "w") as out_file:
        out_file.write(d)


def rmd_create_fig(fname, rname, **kwargs):
    kwargs.setdefault("echo", "FALSE")
    kwargs.setdefault("message", "FALSE")
    kwargs.setdefault("fig.aling", "'center'")
    kwargs.setdefault("fig.caption", "''")
    kwargs.setdefault("out.width", "'0.75\\linewidth'")
    kwargs.setdefault("fig.pos", "'H'")

    opts = [f"{key}={value}" for key, value in kwargs.items()]

    return """
    ```r {rname}, ",".join(opts)
    knitr::include_graphics({fname}")
    ```
    """


def report_read_summary_template_input(wildcards):
  targets = {
      "rmd": f"{workflow.basedir}/report/read_summary.Rmd",
      "regions": REGIONS,
      "plot_total": join_path("results/plots/read_summary/total.pdf")}
  with open(REGIONS, "r") as f:
    regions = f.readlines()

  targets += {"plot_{region}": join_path("results/plots/read_summary", f"{region}.pdf")
             for region in regions}

  return targets


rule report_read_summary_template:
  input: report_read_summary_template_input
  params:
    mapping={"OUTPUT": join_path("/")},
  run:
    mapping = params.mapping + rmd_add_plot(input.plot_total)

    with open(input["rmd"], "r") as in_file:
      s = in_file.read()
      d = template(s, params.mapping)
      with open(output[0], "w") as out_file:
        out_file.write(d)


rule report_feature_summary_template:
  input: f"{workflow.basedir}/report/feature_summary.Rmd",
  output: join_path("results/report/feature_summary.Rmd")
  params:
    mapping={"OUTPUT": join_path("results/report/")},
  run:
    with open(input[0], "r") as in_file:
      s = in_file.read()
      d = template(s, params.mapping)
      with open(output[0], "w") as out_file:
        out_file.write(d)


rule report_session_template:
  input: f"{workflow.basedir}/report/session.Rmd",
  output: join_path("results/report/session.Rmd")
  params:
    mapping={"OUTPUT": join_path("results/report/")},
  run:
    with open(input[0], "r") as in_file:
      s = in_file.read()
      d = template(s, params.mapping)
      with open(output[0], "w") as out_file:
        out_file.write(d)
