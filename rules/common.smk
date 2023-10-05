import os
import re
from snakemake.shell import shell

wildcard_constraints:
  ANALYSIS = "(analysis|downsampling)",


ruleorder: plot_stats_sn_summary > plot_stats_section_summary


SAMPLES = pep.sample_table
PRONTO = pep.config.pronto


def join_path(*e):
  return os.path.join(PRONTO.get("output_dir"), *e)


REF_FASTA = join_path("data/ref.fasta")
MODS = join_path("data/mods.tsv")
REGIONS = join_path("results/data/regions.txt")


STATS_SECTION2COLUMNS = {
  "RL": ["read_length", "count", ],
  "MAPQ": ["mapq", "count", ],
  "ID": ["length", "insertion", "deletion", ],
}


STATS_SN_COLUMNS = ["reads mapped",
          "reads MQ0",
          "reads QC failed",
          "error rate",
          "average length",
          "average quality",]
STATS_SN_COLUMN2LABEL = {col.replace(" ", "_"): col for col in STATS_SN_COLUMNS}


DIR_READ_SUMMARY_PLOTS = join_path("plots/read_summary")
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


# TODO move to common
def read_summary_rdfs(wildcards):
  fnames = FNAMES_READ_SUMMARY_PLOTS
  fnames += [join_path(f"plots/samtools/stats/merged_SN_{col.replace(' ', '_')}_summary.rds")
             for col in STATS_SN_COLUMNS]
  fnames += [join_path(f"plots/samtools/stats/merged_{col}_summary.rds")
             for col in ["RL", "MAPQ", "I", "D"]]

  return fnames


def analysis_targets():
  targets = analysis_lof_results()
  targets.extend(analysis_feature_plots())

  return targets


def downsampling_targets():
  if not "downsampling" in config:
    return ""

  targets = downsampling_lof_results()

  return targets


def auto_targets():
  key2callback = {
    "downsampling": downsampling_targets,
  }
  targets = analysis_targets()
  for key, callback in key2callback.items():
    if key in config:
      targets.extend(callback())

  targets.append(join_path("results/merged_lof.tsv"))
  targets.append(join_path("results/samtools/stats/merged_SN.tsv"))
  # FIXME targets.append(join_path("results/plots/read_length_summary.pdf"))

  targets.append(join_path("plots/mod_summary.pdf"))
  targets.append(join_path("plots/read_summary/total.pdf"))
  targets.append(join_path("plots/feature_lof_summary.pdf"))
  targets.append(join_path("plots/analysis/feature_summary.pdf"))
  targets.append(join_path("report/report.html"))

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
          join_path("plots/feature/analysis",
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

  return targets


rule merge_lof_results:
  input: merged_lof_results()
  output: join_path("results/merged_lof.tsv"),
  log: join_path("logs/merge_lof_results.log"),
  run:
    def helper(fname):
      df = pd.read_csv(fname, sep="\t")
      df["fname"] = fname
      result = re.search(r".+/results/([^/]+)/jacusa2/([^/]+)/lof/neighbors~([0-9]+)_contamination~([0-9.,]+)/[^/]+", fname)
      analysis, bam_prefix, neighbors, contamination = result.groups()
      df["analysis"] = analysis
      df["parameters"] = bam_prefix
      df["lof_neightbors"] = neighbors
      df["lof_contamination"] = contamination
      if "seed~" in bam_prefix and "reads~" in bam_prefix:
        result = re.search(r"seed~(.+)_reads~([0-9]+)", bam_prefix)
        seed, reads = result.groups()
        df["downsampling_seed"] = seed
        df["downsampling_reads"] = reads
      return df

    dfs = [helper(fname) for fname in input]
    df = pd.concat(dfs, ignore_index=True)
    df.to_csv(output[0], sep="\t", index=False)
