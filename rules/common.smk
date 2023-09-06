import os

from snakemake.shell import shell

wildcard_constraints:
  ANALYSIS = "(analysis|downsampling|mixing)"


SAMPLES = pep.sample_table
PRONTO = pep.config.pronto


def join_path(*e):
  return os.path.join(PRONTO.get("output_dir"), *e)


REF_FASTA = join_path("data/ref.fasta")
MODS = join_path("data/mods.tsv")
REGIONS = join_path("results/data/regions.txt")


def analysis_targets():
  targets = []
  for lof in config["lof"]:
    neighbors = lof["neighbors"]
    contamination = lof["contamination"]
    targets.append(
        join_path("results/analysis/jacusa2/preprocessed/lof",
                  f"neighbors~{neighbors}_contamination~{contamination}",
                  "cond1_vs_cond2_lof.tsv"))

  print(targets)
  return targets



def downsampling_targets():
  if not "downsampling" in config:
    return ""

  targets = []
  path = "resuts/downsampling/bams"
  for seed in config["downsampling"]["seeds"]:
    for reads in config.downsampling.reads:
      targets.append(join_path(path,
                               f"seed~{seed}_reads~{reads}",
                               "cond1_vs_cond2"))

  return [targets]


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


def create_include_bam_rules(condition):
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


create_include_bam_rules(1)
create_include_bam_rules(2)


rule create_regions:
  output: REGIONS
  params:
    regions=PRONTO.regions
  run:
      with open(output[0], "w") as f:
        for region in set(params.regions):
          f.write(f"{region}\n")


def replicates(condition):
  return len(SAMPLES.loc[SAMPLES["condition"] == str(condition), "filename"].tolist())


rule link_analysis_bams:
  input: join_path("results/data/bams/preprocessed/{filename}.bam"),
  output: join_path("results/analysis/bams/preprocessed/{filename}.bam"),
  shell: """
    ln -sf ../../../data/bams/preprocessed/{wildcards.filename}.bam {output}
  """
