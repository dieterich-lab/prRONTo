import jinja2

# FIXME remove
def template(s, mapping):
  return s


rule report_render:
  input: report=join_path("results/report/report.Rmd"),
         config=join_path("results/report/config.Rmd"),
         read_summary=join_path("results/report/read_summary.Rmd"),
         feature_summary=join_path("results/report/feature_summary.Rmd"),
         session=join_path("results/report/session.Rmd"),
         params=join_path("results/report/params.yaml"),
  output: join_path("results/report/report.html"),
  log: join_path("logs/report/render.log"),
  shell: """
    Rscript {workflow.basedir}/scripts/render_report.R \
      --params {input.params} \
      --format html_document \
      --output {output} \
      {input.report} \
      2> {log}
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


def parse_template(input, params, output):
  mapping = dict(params)
  with open(input["rmd"], "r") as in_file:
    s = in_file.read()
    t = jinja2.Template(s)
    r = t.render(**mapping)
    with open(output[0], "w") as out_file:
      out_file.write(r)


rule report_read_summary_template:
  input: rmd=f"{workflow.basedir}/report/read_summary.Rmd",
         regions=REGIONS,
         rds=FNAMES_READ_SUMMARY_PLOTS +
             [join_path("results/plots/read_length_summary.rds"),
              join_path("results/plots/read_mapq_summary.rds"),],
  output: join_path("results/report/read_summary.Rmd")
  params:
    config=dict(config),
    PRONTO=dict(PRONTO),
  run:
    parse_template(input, params, output)


rule report_feature_summary_template:
  input: rmd=f"{workflow.basedir}/report/feature_summary.Rmd",
         rds=[fname for fname in FNAMES_READ_SUMMARY_PLOTS if fname.endswith("rds")]
  output: join_path("results/report/feature_summary.Rmd")
  run:
    with open(input["rmd"], "r") as in_file:
      s = in_file.read()
      d = template(s, params)
      with open(output[0], "w") as out_file:
        out_file.write(d)


rule report_create_params:
  output: join_path("results/report/params.yaml")
  params:
    config=config,
    pronto=dict(PRONTO),
  run:
    d = dict(params)
    d["meta"] = {
      "configfile": workflow.configfiles,
      "pepfile": config["pepfile"],
      "workdir": os.getcwd(),
    }
    with open(output[0], "w") as f:
       yaml.safe_dump(d, f)


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
