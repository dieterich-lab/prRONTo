rule report_create:
  input: f"{workflow.basedir}/scripts/create_report.Rmd"
  output: join_path("results/report/report.pdf"),
  log: join_path("logs/report/create.log"),
  shell: """
    Rscript {workflow.basedir}/scripts/create_report.R --output {output} {input}
  """
