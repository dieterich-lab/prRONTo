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
  log: join_path("logs/{ANALYSIS}/add_meta/{bam_prefix}/{comparison}.log"),
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


def analysis_lof_results():
  targets = []
  for lof in config["lof"]:
    neighbors = lof["neighbors"]
    contamination = lof["contamination"]
    targets.append(
        join_path("results/analysis/jacusa2",
                  "preprocessed",
                  "lof",
                  f"neighbors~{neighbors}_contamination~{contamination}",
                  "cond1_vs_cond2.tsv"))

  return targets


def downsampling_lof_results():
  targets = []
  for seed in config["downsampling"]["seed"]:
    for lof in config["lof"]:
      neighbors = lof["neighbors"]
      contamination = lof["contamination"]
      for reads in config["downsampling"]["reads"]:
        targets.append(
          join_path("results/downsampling/jacusa2",
                    f"seed~{seed}_reads~{reads}",
                    "lof",
                    f"neighbors~{neighbors}_contamination~{contamination}",
                    "cond1_vs_cond2.tsv"))

  return targets


def merged_lof_results():
  targets = analysis_lof_results()

  if "downsampling" in config:
    targets.extend(downsampling_lof_results())

    if "mixing" in config:
      targets.extend(downsampling_lof_results())

  return targets


rule merge_lof_results:
  input: merged_lof_results()
  output: join_path("results/merged_lof.tsv"),
  run:
    dfs = [pd.read_csv(fname, sep="\t") for fname in input]
    df = pd.concat(dfs, ignore_index=True, )
    df.to_csv(output[0], sep="\t", index=False)
