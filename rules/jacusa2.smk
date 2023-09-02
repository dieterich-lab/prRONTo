def replicates(cond):
  # TODO do calculation
  return 1

  if wildcards.analysis == "analysis":
    analysis_prefix = "results/data/preprocessed_bams"
  else:
    analysis_prefix = f"results/{wildcards.analysis}/bams/bams/bams/bams/bams/bams/bams"
   cond{cond}_rep{rep}.bam/


def bams(wildcards, prefix):
  bams = []
  for rep in xrange(replicates(cond)):
    bams.append(join_path(""))

  return bams


def bais(cond):
  return [f"{bam}.bai" for bam in bams_analysis(cond)]


rule jacusa2_call_analysis:
  input:
    bamsA=bams_analysis(1),
    bais1=bais_analysis(1),
    bamsA=bams_analysis(2),
    bais2=bais_analysis(2),
  output: join_path("results/analysis/jacusa2/{condA}_vs_{condB}.out")
  log: join_path("logs/analysis/jacusa2/preprocessed_bams/{condA}_vs_{condB}.log")
  threads: config["jacusa2"]["threads"]
  params:
    opts=config["jacusa2"]["opts"],
    min_mapq=config["jacusa2"]["min_mapq"],
    min_bq=config["jacusa2"]["min_bq"],
    min_cov=config["jacusa2"]["min_cov"],
    input=lambda wildcards, input: ",".join(input.bamsA) + " " + ",".join(input.bamsB),
  shell: """
    jacusa2 call-2 \
        {params.opts} \
        -m {params.min_mapq} \
        -q {params.min_bq} \
        -c {params.min_cov} \
        -p {threads} \
        -r {output} \
        {params.input} \
        2> {log}
  """


rule jacusa2_add_scores_analysis:
  input: join_path("results/analysis/jacusa2/{condA}_vs_{condB}.out")
  output: join_path("results/analysis/jacusa2/{condA}_vs_{condB}_scores.tsv")
  log: join_path("results/jacusa2/add_scores/{ANALYSIS}/{bams}/{condA}_vs_{condB}.log")
  shell: """
    Rscript --vanilla {workflow.basedir}/scripts/add_scores.R \
        --output {output} \
        --comparison {wildcards.condA}_vs_{wildcards.condB} \
        --analysis {wildcards.ANALYSIS} \
        --parameters {wildcards.bams} \
        {input} \
        2> {log}
  """

