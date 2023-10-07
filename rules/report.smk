import jinja2 # TODO use snakemake 7 engine


rule report_render:
  input: report=join_path("report/report.Rmd"),
         config=join_path("report/config.Rmd"),
         read_summary=join_path("report/read_summary.Rmd"),
         feature_summary=join_path("report/feature_summary.Rmd"),
         session=join_path("report/session.Rmd"),
         params=join_path("report/params.yaml"),
  output: join_path("report/report.html"),
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
  input: rmd=f"{workflow.basedir}/report/report.Rmd",
  output: join_path("report/report.Rmd")
  params:
    config=dict(config),
    PRONTO=dict(PRONTO),
  run:
    parse_template(input, params, output)


rule report_config_template:
  input: rmd=f"{workflow.basedir}/report/config.Rmd",
         rds=[join_path("plots/mod_summary.rds"), ],
  output: join_path("report/config.Rmd")
  params:
    config=dict(config),
    PRONTO=dict(PRONTO),
  run:
    parse_template(input, params, output)


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
         rds=read_summary_rdfs,
  output: join_path("report/read_summary.Rmd")
  params:
    config=dict(config),
    PRONTO=dict(PRONTO),
  run:
    parse_template(input, params, output)


# FIXME rds
rule report_feature_summary_template:
  input: rmd=f"{workflow.basedir}/report/feature_summary.Rmd",
         rds=[fname for fname in FNAMES_READ_SUMMARY_PLOTS if fname.endswith("rds")] +
             [join_path("plots/analysis/feature_summary.rds"), ]
  output: join_path("report/feature_summary.Rmd")
  params:
    config=dict(config),
    PRONTO=dict(PRONTO),
  run:
    parse_template(input, params, output)


rule report_create_params:
  output: join_path("report/params.yaml")
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
  input: rmd=f"{workflow.basedir}/report/session.Rmd",
  output: join_path("report/session.Rmd")
  params:
    config=dict(config),
    PRONTO=dict(PRONTO),
  run:
    parse_template(input, params, output)
