import jinja2


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
    pep=pep,
    basedir=workflow.basedir
  run:
    parse_template(input, params, output)


rule report_config_template:
  input: rmd=f"{workflow.basedir}/report/config.Rmd",
         rds=[join_path("plots/mod_summary.rds"), ],
  output: join_path("report/config.Rmd")
  params:
    config=dict(config),
    PRONTO=dict(PRONTO),
    basedir=workflow.basedir
  run:
    parse_template(input, params, output)


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
    basedir=workflow.basedir
  run:
    parse_template(input, params, output)


# FIXME rds **
# path and expand (ANALYSIS, bam_prefix, neighbors, contamination, comparison, feature, seq_id, suffix
rule report_feature_summary_template:
  input: rmd=f"{workflow.basedir}/report/feature_summary.Rmd",
         rds=[original_feature_plots()] +
             [downsampling_feature_plots()] +
             [downsampling_summary_plots()] +
             [join_path("plots/original/feature_summary.rds"), ]
  output: join_path("report/feature_summary.Rmd")
  params:
    config=dict(config),
    PRONTO=dict(PRONTO),
    basedir=workflow.basedir
  run:
    parse_template(input, params, output)


rule report_session_template:
  input: rmd=f"{workflow.basedir}/report/session.Rmd",
  output: join_path("report/session.Rmd")
  params:
    config=dict(config),
    PRONTO=dict(PRONTO),
    basedir=workflow.basedir
  run:
    parse_template(input, params, output)


rule report_create_params:
  output: join_path("report/params.yaml")
  params:
    config=config,
    PRONTO=dict(PRONTO),
    basedir=workflow.basedir
  run:
    d = dict(params)
    d["meta"] = {
      "configfile": workflow.configfiles,
      "pepfile": config["pepfile"],
      "workdir": os.getcwd(),
    }
    with open(output[0], "w") as f:
       yaml.safe_dump(d, f)
