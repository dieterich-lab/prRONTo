def jacusa2_bams(condition):
  def helper(wildcards):
    condition2wildcard = {
        "A": wildcards.condA,
        "B": wildcards.condB,}
    bams = []
    bams_pattern = join_path("results",
                             wildcards.ANALYSIS,
                             "bams",
                             wildcards.bam_prefix,
                             condition2wildcard[condition] + "_rep{rep}.bam")
    # FIXME change to condI to enable multiple combinations
    if condition2wildcard[condition] == "cond1":
      reps = REPLICATES1
    elif condition2wildcard[condition] == "cond2":
      reps = REPLICATES2
    else:
      reps = 1

    for rep in range(reps):
      bams.append(bams_pattern.format(rep = str(rep + 1)))

    return bams

  return helper


def jacusa2_bais(condition):
  def helper(wildcards):
    return [f"{bam}.bai" for bam in jacusa2_bams(condition)(wildcards)]

  return helper


# FIXME make generic comparisons
# FIXME make library type a config option
rule jacusa2_call:
  input: # condA can be cond1,2, or 3 from the config - therefore using A, B
    bamsA=jacusa2_bams("A"),
    baisA=jacusa2_bais("A"),
    bamsB=jacusa2_bams("B"),
    baisB=jacusa2_bais("B"),
  output: join_path("results/{ANALYSIS}/jacusa2/{bam_prefix}/{condA}_vs_{condB}.out"),
  log: join_path("logs/{ANALYSIS}/jacusa2/call2/{bam_prefix}/{condA}_vs_{condB}.log"),
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


def jacusa2_features():
  return ",".join(config["jacusa2"]["features"])


rule jacusa2_add_scores:
  input: join_path("results/{ANALYSIS}/jacusa2/{bam_prefix}/{comparison}.out"),
  output: temp(join_path("results/{ANALYSIS}/jacusa2/{bam_prefix}/{comparison}_scores.tsv")),
  log: join_path("logs/{ANALYSIS}/jacusa2/add_scores/{bam_prefix}/{comparison}.log"),
  params:
    features=jacusa2_features(),
  shell: """
    Rscript --vanilla {workflow.basedir}/scripts/add_features.R \
        --output {output} \
        --features {params.features} \
        --comparison {wildcards.comparison} \
        --analysis {wildcards.ANALYSIS} \
        --parameters {wildcards.bam_prefix} \
        {input} \
        2> {log}
  """
