rule plot_mod_summary:
  input: join_path("data/mods.tsv"),
  output: plot=join_path("plots/mod_summary.pdf"),
          rds=join_path("plots/mod_summary.rds"),
  log: join_path("logs/plot/mod_summary.log"),
  shell: """
    Rscript {workflow.basedir}/scripts/plot_mod_summary.R \
      --output {output.plot} \
      {input} \
      2> {log}
  """


# TODO replicate violin plot with table
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
  output: plot=join_path("plots/feature_lof_summary.pdf"),
          rdf=join_path("plots/feature_lof_summary.rds"),
  log: join_path("logs/plot/feature_lof_summary.log"),
  shell: """
    Rscript {workflow.basedir}/scripts/plot_feature_lof_summary.R \
      --output {output.plot} \
      {input} \
      2> {log}
  """


# TODO log path add feature
rule plot_feature:
  input: join_path("results/{ANALYSIS}/jacusa2/{bam_prefix}/lof/neighbors~{neighbors}_contamination~{contamination}/{comparison}.tsv"),
  output: directory(join_path("plots/feature/{ANALYSIS}/{bam_prefix}/neighbors~{neighbors}_contamination~{contamination}/{comparison}/feature~{feature}")),
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


rule plot_stats_section_summary:
  input: join_path("results/samtools/stats/merged_{section}.tsv"),
  output: pdf=join_path("plots/samtools/stats/merged_{section}_summary.pdf"),
          rds=join_path("plots/samtools/stats/merged_{section}_summary.rds"),
  log: join_path("logs/plot/samtools/stats/merged_{section}_summary.log"),
  params:
    xvalue=lambda wildcards: STATS_SECTION2COLUMNS[wildcards.section][0],
    yvalue=lambda wildcards: STATS_SECTION2COLUMNS[wildcards.section][1]
  shell: """
    Rscript {workflow.basedir}/scripts/plot_stats_section_summary.R \
      --output {output.pdf} \
      --xvalue "{params.xvalue}" \
      --yvalue "{params.yvalue}" \
      {input} \
      2> {log}
  """


use rule plot_stats_section_summary as plot_stats_I_summary with:
  input: join_path("results/samtools/stats/merged_ID.tsv"),
  output: pdf=join_path("plots/samtools/stats/merged_I_summary.pdf"),
          rds=join_path("plots/samtools/stats/merged_I_summary.rds"),
  log: join_path("logs/plot/samtools/stats/merged_I_summary.log"),
  params:
    xvalue=lambda wildcards: STATS_SECTION2COLUMNS["ID"][0],
    yvalue=lambda wildcards: STATS_SECTION2COLUMNS["ID"][1]


use rule plot_stats_section_summary as plot_stats_D_summary with:
  input: join_path("results/samtools/stats/merged_ID.tsv"),
  output: pdf=join_path("plots/samtools/stats/merged_D_summary.pdf"),
          rds=join_path("plots/samtools/stats/merged_D_summary.rds"),
  log: join_path("logs/plot/samtools/stats/merged_D_summary.log"),
  params:
    xvalue=lambda wildcards: STATS_SECTION2COLUMNS["ID"][0],
    yvalue=lambda wildcards: STATS_SECTION2COLUMNS["ID"][2]


rule plot_stats_sn_summary:
  input: join_path("results/samtools/stats/merged_SN.tsv"),
  output: pdf=join_path("plots/samtools/stats/merged_SN_{column}_summary.pdf"),
          rds=join_path("plots/samtools/stats/merged_SN_{column}_summary.rds"),
  log: join_path("logs/plot/samtools/stats/merged_SN_{column}_summary.log"),
  params:
    label=lambda wildcards: STATS_SN_COLUMN2LABEL[wildcards.column]
  shell: """
    Rscript {workflow.basedir}/scripts/plot_stats_sn_summary.R \
      --output {output.pdf} \
      --column "{params.label}" \
      --label "{params.label}" \
      {input} \
      2> {log}
  """


rule plot_feature_summary:
  input: join_path("results/merged_lof.tsv"),
  output: pdf=join_path("plots/{ANALYSIS}/feature_summary.pdf"),
          rds=join_path("plots/{ANALYSIS}/feature_summary.rds"),
  log: join_path("logs/plot/{ANALYSIS}/feature_summary.log"),
  shell: """
    Rscript {workflow.basedir}/scripts/plot_feature_summary.R \
      --output {output.pdf} \
      {input} \
      2> {log}
  """
