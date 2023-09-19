import os
import re
from snakemake.shell import shell

wildcard_constraints:
  ANALYSIS = "(analysis|downsampling|mixing)",


SAMPLES = pep.sample_table
PRONTO = pep.config.pronto


def join_path(*e):
  return os.path.join(PRONTO.get("output_dir"), *e)


REF_FASTA = join_path("data/ref.fasta")
MODS = join_path("data/mods.tsv")
REGIONS = join_path("results/data/regions.txt")


def analysis_targets():
  targets = analysis_lof_results()
  targets.extend(analysis_feature_plots())

  return targets



def downsampling_targets():
  if not "downsampling" in config:
    return ""

  targets = downsampling_lof_results()

  return targets


def mixing_targets():
  if not "mixing" in config:
    return ""

  return [] # pass


def auto_targets():
  key2callback = {
    "mixing": mixing_targets,
    "downsampling": downsampling_targets,
  }
  targets = analysis_targets()
  for key, callback in key2callback.items():
    if key in config:
      targets.extend(callback())

  targets.append(join_path("results/merged_lof.tsv"))

  targets.append(join_path("results/plots/mod_summary.pdf"))
  targets.append(join_path("results/plots/read_summary/total.pdf"))
  targets.append(join_path("results/plots/feature_lof_summary.pdf"))
  targets.append(join_path("results/report/report.html"))

  return targets


rule include_fasta:
  input: PRONTO.ref
  output: REF_FASTA
  params:
    include=config.get("include", {}).get("ref")
  run:
      if params.include == "copy":
       cmd = "cp"
      else:
       cmd = "ln -s"

      shell(cmd + " {input} {output}")


rule include_mods:
  input: PRONTO.mods
  output: MODS
  params:
    include=config.get("include", {}).get("mods")
  run:
      if params.include == "copy":
        cmd = "cp"
      else:
        cmd = "ln -s"

      shell(cmd + " {input} {output}")


# FIXME what if there is a replicate column in SAMPLES
def create_include_bam_rules(condition: int):
  fnames = SAMPLES.loc[SAMPLES["condition"] == str(condition), "filename"]
  for i, fname in enumerate(fnames, start=1):
    rule:
      name: f"dyn_include_bam_cond{condition}_rep{i}"
      input: fname
      output: join_path(f"data/bams/raw/cond{condition}_rep{i}.bam")
      params:
        include=config.get("include", {}).get("bams")
      run:
        if params.include == "copy":
          cmd = "cp"
        else:
          cmd = "ln -s"

        shell(cmd + " {input} {output}")


# create dynamically rules to include raw bams for all conditions
for condition in SAMPLES["condition"].unique():
  create_include_bam_rules(condition)


rule create_regions:
  output: REGIONS
  params:
    regions=PRONTO.regions
  run:
      with open(output[0], "w") as f:
        for region in set(params.regions):
          f.write(f"{region}\n")


def get_replicates(condition):
  return len(SAMPLES.loc[SAMPLES["condition"] == str(condition), "filename"].tolist())


rule link_analysis_bams:
  input: join_path("results/data/bams/preprocessed/{filename}.bam"),
  output: join_path("results/analysis/bams/preprocessed/{filename}.bam"),
  shell: """
    ln -sf ../../../data/bams/preprocessed/{wildcards.filename}.bam {output}
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

# TODO add region
def analysis_feature_plots():
  targets = []
  for lof in config["lof"]:
    neighbors = lof["neighbors"]
    contamination = lof["contamination"]
    for feature in config["jacusa2"]["features"]:
      targets.append(
          join_path("results/plots/feature/analysis",
                    "preprocessed",
                    f"neighbors~{neighbors}_contamination~{contamination}",
                    "cond1_vs_cond2",
                    f"feature~{feature}"))

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
      pass # TODO mixing

  return targets


rule merge_lof_results:
  input: merged_lof_results()
  output: join_path("results/merged_lof.tsv"),
  log: join_path("logs/merge_lof_results.log"), # TODO - where to put this
  run:
    dfs = [pd.read_csv(fname, sep="\t") for fname in input]
    df = pd.concat(dfs, ignore_index=True, )
    df.to_csv(output[0], sep="\t", index=False)
