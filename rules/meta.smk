import re


def get_max_context():
  features = config["jacusa2"]["features"]
  contexts = [5, ]
  for feature in features:
    contexts.extend(re.findall("\d+", feature))
  contexts = [int(context) for context in contexts]

  return max(contexts)


rule meta_add:
  input: scores=join_path("results/{ANALYSIS}/jacusa2/{bam_prefix}/{comparison}_scores.tsv"),
         ref=REF_FASTA,
         mods=MODS,
  output: temp(join_path("results/{ANALYSIS}/jacusa2/{bam_prefix}/{comparison}_meta.tsv")),
  log: join_path("logs/{ANALYSIS}/add_meta/{bam_prefix}/{comparison}.log"),
  params: context=get_max_context()
  shell: """
    Rscript {workflow.basedir}/scripts/add_meta.R \
        --output {output} \
        --mods {input.mods} \
        --fasta {input.ref} \
        --context {params.context} \
        {input.scores} \
        2> {log}
  """


def parse_features():
  opts = [f"--feature {feature}" for feature in config["jacusa2"]["features"]]

  return " ".join(opts)


rule meta_add_lof:
  input: join_path("results/{ANALYSIS}/jacusa2/{bam_prefix}/{comparison}_meta.tsv"),
  output: join_path("results/{ANALYSIS}/jacusa2/{bam_prefix}/lof/neighbors~{neighbors}_contamination~{contamination}/{comparison}.tsv"),
  log: join_path("logs/{ANALYSIS}/add_lof/{bam_prefix}/neighbors~{neighbors}_contamination~{contamination}/{comparison}.log"),
  params:
    features=parse_features()
  shell: """
    python {workflow.basedir}/scripts/add_lof.py \
        {params.features} \
        --filter-median \
        --neighbors {wildcards.neighbors} \
        --contamination {wildcards.contamination} \
        --output {output} \
        {input} \
        2> {log}
  """
