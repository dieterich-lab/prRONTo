def _plot_mods_opts(wildcards):
  if wildcards.region == "summary":
    return ""

  return f"-s {wildcards.region}"

rule plot_mod_summary:
  input: join_path("data/mods.tsv"),
  output: pdf=join_path("plots/mods/{region}.pdf"),
          rds=join_path("plots/mods/{region}.rds"),
  log: join_path("logs/plot/mods/{region}.log"),
  params:
    opts=_plot_mods_opts,
  shell: """
    Rscript {workflow.basedir}/scripts/plot_mod_summary.R \
      --output {output.pdf} \
      {params.opts} \
      {input} \
      2> {log}
  """


rule plot_read_summary:
  input: join_path("results/coverage_summary.tsv"),
  output: files=FNAMES_READ_SUMMARY_PLOTS,
          dir=directory(DIR_READ_SUMMARY_PLOTS),
  log: join_path("logs/plot/read_summary.log"),
  params:
    cond1=PRONTO["condition1"],
    cond2=PRONTO["condition2"],
  shell: """
    Rscript {workflow.basedir}/scripts/plot_read_summary.R \
      --cond1 {params.cond1} \
      --cond2 {params.cond2} \
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


def dir_plot_feature():
  return join_path("plots/feature/{ANALYSIS}/{bam_prefix}/neighbors~{neighbors}_contamination~{contamination}/feature~{feature}")


def files_plot_feature(suffix):
  return [dir_plot_feature() + f"/{seq_id}{suffix}" for seq_id in PRONTO["regions"]]


rule plot_feature:
  input: join_path("results/{ANALYSIS}/jacusa2/{bam_prefix}/lof/neighbors~{neighbors}_contamination~{contamination}/cond1_vs_cond2.tsv"),
  output: dir=directory(dir_plot_feature()),
          pdf=files_plot_feature(".pdf"),
          rds=files_plot_feature(".rds") + files_plot_feature("_barplot.rds"),
          tsv=files_plot_feature("_outlier.tsv"),
  log: join_path("logs/plot/{ANALYSIS}/{bam_prefix}/neighbors~{neighbors}_contamination~{contamination}/feature~{feature}.log"),
  params:
    targets="", # TODO add targets
  shell: """
    Rscript {workflow.basedir}/scripts/plot_feature.R \
      {params.targets} \
      --output {output.dir} \
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
    yvalue=lambda wildcards: STATS_SECTION2COLUMNS[wildcards.section][1],
    xlabel=lambda wildcards: STATS_SECTION2COLUMNS[wildcards.section][0].replace("_", " "),
    cond1=PRONTO["condition1"],
    cond2=PRONTO["condition2"],
  shell: """
    Rscript {workflow.basedir}/scripts/plot_stats_section_summary.R \
      --output {output.pdf} \
      --xvalue "{params.xvalue}" \
      --yvalue "{params.yvalue}" \
      --xlabel "{params.xlabel}" \
      --cond1 "{params.cond1}" \
      --cond2 "{params.cond2}" \
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
    yvalue=lambda wildcards: STATS_SECTION2COLUMNS["ID"][1],
    xlabel=lambda wildcards: STATS_SECTION2COLUMNS["ID"][0].replace("_", " "),
    cond1=PRONTO["condition1"],
    cond2=PRONTO["condition2"],


use rule plot_stats_section_summary as plot_stats_D_summary with:
  input: join_path("results/samtools/stats/merged_ID.tsv"),
  output: pdf=join_path("plots/samtools/stats/merged_D_summary.pdf"),
          rds=join_path("plots/samtools/stats/merged_D_summary.rds"),
  log: join_path("logs/plot/samtools/stats/merged_D_summary.log"),
  params:
    xvalue=lambda wildcards: STATS_SECTION2COLUMNS["ID"][0],
    yvalue=lambda wildcards: STATS_SECTION2COLUMNS["ID"][2],
    xlabel=lambda wildcards: STATS_SECTION2COLUMNS["ID"][0].replace("_", " "),
    cond1=PRONTO["condition1"],
    cond2=PRONTO["condition2"],


def _norm(wildcards):
  if wildcards.column == "reads_mapped":
    return '--normalize "raw total sequences"'

  return ""


rule plot_stats_sn_summary:
  input: join_path("results/samtools/stats/merged_SN.tsv"),
  output: pdf=join_path("plots/samtools/stats/merged_SN_{column}_summary.pdf"),
          rds=join_path("plots/samtools/stats/merged_SN_{column}_summary.rds"),
  log: join_path("logs/plot/samtools/stats/merged_SN_{column}_summary.log"),
  params:
    label=lambda wildcards: STATS_SN_COLUMN2LABEL[wildcards.column],
    cond1=PRONTO["condition1"],
    cond2=PRONTO["condition2"],
    norm=_norm,
  shell: """
    Rscript {workflow.basedir}/scripts/plot_stats_sn_summary.R \
      {params.norm} \
      --cond1 {params.cond1} \
      --cond2 {params.cond2} \
      --output {output.pdf} \
      --column "{params.label}" \
      --label "{params.label}" \
      {input} \
      2> {log}
  """


rule plot_feature_summary:
  input: join_path("results/merged_lof.tsv"),
  output: pdf=join_path("plots/{ANALYSIS}/feature_summary/neighbors~{neighbors}_contamination~{contamination}.pdf"),
          rds=join_path("plots/{ANALYSIS}/feature_summary/neighbors~{neighbors}_contamination~{contamination}.rds"),
  log: join_path("logs/plot/{ANALYSIS}/feature_summary/neighbors~{neighbors}_contamination~{contamination}.log"),
  shell: """
    Rscript {workflow.basedir}/scripts/plot_feature_summary.R \
      --lof_params {wildcards.neighbors}:{wildcards.contamination} \
      --output {output.pdf} \
      {input} \
      2> {log}
  """


rule plot_downsampling_summary:
  input: join_path("results/merged_lof.tsv"),
  output: pdf=join_path("plots/downsampling_summary/neighbors~{neighbors}_contamination~{contamination}/feature~{feature}/{region}.pdf")
  log: join_path("logs/plot/downsampling_summary/neighbors~{neighbors}_contamination~{contamination}/feature~{feature}/{region}.log"),
  params:
    lof_params=lambda wildcards: f"{wildcards.neighbors}:{wildcards.contamination}"
  shell: """
    Rscript {workflow.basedir}/scripts/plot_upset.R \
      --features {wildcards.feature} \
      --lof_params {params.lof_params} \
      --regions {wildcards.region} \
      --output {output.pdf} \
      {input} \
      2> {log}
  """
