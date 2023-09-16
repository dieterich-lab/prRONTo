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


# TODO plot reads per seq_id contig
# replicate violin plot with table
# FIXME reads 200 -> in plot 400
rule plot_read_summary:
  input: join_path("results/read_summary.tsv"),
  output: total=join_path("results/plots/read_summary/total.pdf"),
          dir=directory(join_path("results/plots/read_summary")),
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
