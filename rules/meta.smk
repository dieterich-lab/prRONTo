import re


def get_max_context():
  features = config["jacusa2"]["features"]
  contexts = [5, ]
  for feature in features:
    context = re.sub(r"[MDI_]+", "", feature)
    if context:
      contexts.append(int(context))

  return max(contexts)


rule add_meta:
  input: scores=join_path("results/{ANALYSIS}/jacusa2/{bam_prefix}/{comparison}_scores.tsv"),
         ref=REF_FASTA,
         mods=MODS,
  output: temp(join_path("results/{ANALYSIS}/jacusa2/{bam_prefix}/{comparison}_meta.tsv")),
  log: join_path("logs/add_meta/{ANALYSIS}/{bam_prefix}/{comparison}.log"),
  params:
    context=get_max_context()
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


rule add_lof:
  input: join_path("results/{ANALYSIS}/jacusa2/{bam_prefix}/{comparison}_meta.tsv"),
  output: join_path("results/{ANALYSIS}/jacusa2/{bam_prefix}/lof/neighbors~{neighbors}_contamination~{contamination}/{comparison}_lof.tsv"),
  log: join_path("logs/add_lof/{ANALYSIS}/{bam_prefix}/neighbors~{neighbors}_contamination~{contamination}/{comparison}.log"),
  params:
    features=parse_features()
  shell: """
    python {workflow.basedir}/scripts/add_lof.py \
        --feature {params.features} \
        --filter_median \
        --neighbors {wildcards.neighbors} \
        --contamination {wildcards.contamination} \
        --output {output} \
        {input} \
        2> {log}
  """
