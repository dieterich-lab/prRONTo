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


def report_template_run(input, output):
  with open(input["yaml"], "r") as f:
    params = yaml.safe_load(f)
    parse_template(input["rmd"], params, output)


rule report_template:
  input: rmd=f"{workflow.basedir}/report/{{template}}.Rmd",
         yaml=join_path("report/params.yaml"),
  output: join_path("report/{template}.Rmd")
  run: report_template_run(input, output)


def mods_rdfs():
  fnames = [join_path("plots/mods/summary.rds"), ]
  for region in PRONTO["regions"]:
    fnames.append(join_path(f"plots/mods/{region}.rds"))

  return fnames


rule report_config_template:
  input: rmd=f"{workflow.basedir}/report/config.Rmd",
         yaml=join_path("report/params.yaml"),
         rds=mods_rdfs(),
  output: join_path("report/config.Rmd")
  run: report_template_run(input, output)


def read_yaml(fname):
  with open(fname, "r") as f:
    return yaml.safe_load(f)


def parse_template(input, params, output):
  mapping = dict(params)

  with open(input, "r") as in_file:
    s = in_file.read()
    jinja2.filters.FILTERS.update(
        {
          "read_yaml": read_yaml,
          "to_yaml": yaml.dump,
         })
    t = jinja2.Template(s)
    r = t.render(**mapping)
    with open(output[0], "w") as out_file:
      out_file.write(r)


rule report_read_summary_template:
  input: rmd=f"{workflow.basedir}/report/read_summary.Rmd",
         yaml=join_path("report/params.yaml"),
         regions=REGIONS,
         rds=read_summary_rdfs,
  output: join_path("report/read_summary.Rmd")
  run: report_template_run(input, output)


def feature_summary_plots():
  fnames = []
  for lof in config["lof"]:
    neighbors = lof["neighbors"]
    contamination = lof["contamination"]
    fnames.append(join_path(f"plots/original/feature_summary/neighbors~{neighbors}_contamination~{contamination}.rds"))

  return fnames

# FIXME rds **
# path and expand (ANALYSIS, bam_prefix, neighbors, contamination, comparison, feature, seq_id, suffix
rule report_feature_summary_template:
  input: rmd=f"{workflow.basedir}/report/feature_summary.Rmd",
         yaml=join_path("report/params.yaml"),
         rds=[original_feature_plots()] + [downsampling_feature_plots()] + [downsampling_summary_plots()] + [feature_summary_plots()]
  output: join_path("report/feature_summary.Rmd")
  run: report_template_run(input, output)


rule report_create_params:
  input: md5_pepfile=config["pepfile"] + ".md5",
         md5_ref=PRONTO["ref"] + ".md5",
         md5_mods=PRONTO["mods"] + ".md5",
         md5_configfile=workflow.configfiles[0] + ".md5",
  output: join_path("report/params.yaml")
  params:
    config=config,
    PRONTO=dict(PRONTO),
    basedir=workflow.basedir,
    VERSIONS=VERSIONS,
    pep=pep.to_dict(),
    meta = {
      "configfile": workflow.configfiles,
      "pepfile": config["pepfile"],
      "workdir": os.getcwd(),},
  run:
    d = dict(params)
    d["md5"] = {}
    for key, fname in input.items():
      if key.startswith("md5"):
        with open(fname, "r") as f:
          d["md5"][key.replace("md5_", "")] = f.read().strip()
    with open(output[0], "w") as f:
      yaml.safe_dump(d, f)
