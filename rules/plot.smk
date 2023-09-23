rule plot_mod_summary:
  input: join_path("data/mods.tsv"),
  output: plot=join_path("results/plots/mod_summary.pdf"),
          rds=join_path("results/plots/mod_summary.rds"),
  log: join_path("logs/plot/mod_summary.log"),
  shell: """
    Rscript {workflow.basedir}/scripts/plot_mod_summary.R \
      --output {output.plot} \
      {input} \
      2> {log}
  """


# TODO replicate violin plot with table
DIR_READ_SUMMARY_PLOTS = join_path("results/plots/read_summary")
READ_SUMMARY_SUBPLOTS = ["data", "seq_ids", ]
if "downsampling" in config:
  READ_SUMMARY_SUBPLOTS.append("downsampling")
FNAMES_READ_SUMMARY_PLOTS = []
for region in ["total"] + PRONTO["regions"]:
  for suffix in ["pdf", "rds"]:
    FNAMES_READ_SUMMARY_PLOTS.append(f"{DIR_READ_SUMMARY_PLOTS}/{region}.{suffix}")
    __subplots = list(READ_SUMMARY_SUBPLOTS)
    if region != "total":
      __subplots.remove("seq_ids")
  for subplot in __subplots:
    FNAMES_READ_SUMMARY_PLOTS.append(f"{DIR_READ_SUMMARY_PLOTS}/{region}_{subplot}.rds")
rule plot_read_summary:
  input: join_path("results/read_summary.tsv"),
  output: files=FNAMES_READ_SUMMARY_PLOTS,
          dir=directory(DIR_READ_SUMMARY_PLOTS),
  log: join_path("logs/plot/read_summary.log"),
  shell: """
    Rscript {workflow.basedir}/scripts/plot_read_summary.R \
      --output {output.dir} \
      {input} \
      2> {log}
  """


rule plot_feature_lof_summary:
  input: join_path("results/merged_lof.tsv"),
  output: plot=join_path("results/plots/feature_lof_summary.pdf"),
          rdf=join_path("results/plots/feature_lof_summary.rds"),
  log: join_path("logs/plot/feature_lof_summary.log"),
  shell: """
    Rscript {workflow.basedir}/scripts/plot_feature_lof_summary.R \
      --output {output.plot} \
      {input} \
      2> {log}
  """


rule plot_feature:
  input: join_path("results/{ANALYSIS}/jacusa2/{bam_prefix}/lof/neighbors~{neighbors}_contamination~{contamination}/{comparison}.tsv"),
  output: directory(join_path("results/plots/feature/{ANALYSIS}/{bam_prefix}/neighbors~{neighbors}_contamination~{contamination}/{comparison}/feature~{feature}")),
  log: join_path("logs/plot/{ANALYSIS}/{bam_prefix}/neighbors~{neighbors}_contamination~{contamination}/{comparison}/feature~{feature}.log"),
  params:
    targets="" # TODO add targets
  shell: """
    Rscript {workflow.basedir}/scripts/plot_feature.R \
      {params.targets} \
      --output {output} \
      --feature {wildcards.feature} \
      {input} \
      2> {log}
  """


rule plot_read_length_summary:
  input: join_path("results/merged_read_length.tsv"),
  output: pdf=join_path("results/plots/read_length_summary.pdf"),
          rds=join_path("results/plots/read_length_summary.rds"),
  log: join_path("logs/plot/read_length_summary.log"),
  shell: """
    Rscript {workflow.basedir}/scripts/plot_read_length_summary.R \
      --output {output.pdf} \
      {input} \
      2> {log}
  """


rule plot_read_mapq_summary:
  input: join_path("results/merged_read_mapq.tsv"),
  output: pdf=join_path("results/plots/read_mapq_summary.pdf"),
          rds=join_path("results/plots/read_mapq_summary.rds"),
  log: join_path("logs/plot/read_mapq_summary.log"),
  shell: """
    Rscript {workflow.basedir}/scripts/plot_read_mapq_summary.R \
      --output {output.pdf} \
      {input} \
      2> {log}
  """


rule plot_read_sn_summary:
  input: join_path("results/merged_read_sn.tsv"),
  output: pdf=join_path("results/plots/read_sn_{column}_summary.pdf"),
          rds=join_path("results/plots/read_sn_{column}_summary.rds"),
  log: join_path("logs/plot/read_sn_{column}_summary.log"),
  shell: """
    Rscript {workflow.basedir}/scripts/plot_read_sn_summary.R \
      --output {output.pdf} \
      --column {wildcards.column} \
      --label {wildcards.column} \
      {input} \
      2> {log}
  """


rule plot_feature_summary:
  input: join_path("results/merged_lof.tsv"),
  output: join_path("results/plots/{ANALYSIS}/feature_summary.pdf"),
  log: join_path("logs/plot/{ANALYSIS}/feature_summary.log"),
  shell: """
    Rscript {workflow.basedir}/scripts/plot_feature_summary.R \
      --output {output} \
      {input} \
      2> {log}
  """
