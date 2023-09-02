import os

from snakemake.shell import shell


def join_path(*e):
  return os.path.join(PRONTO.get("output_dir"), *e)


REF_FASTA = join_path("data/ref.fasta")
MODS = join_path("data/mods.tsv")
REGIONS = join_path("results/data/regions.txt")

SAMPLES = pep.sample_table
PRONTO = pep.config.pronto


def analysis_targets():
  print(config)
  print(SAMPLES)
  print(PRONTO)
  return [join_path("jacusa2/analysis/cond1_vs_cond2.out"),]

# TODO
# * add features/scores
# * add meta
# * plots: # number of sites total


def mixing_targets():
  return [] # pass


def downsampling_targets():
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
  fnames = SAMPLES.loc[SAMPLES["condition"] == condition, "filename"]
  for i, fname in enumerate(fnames, start=1):
    rule include_bam:
      name: f"{dyn_include_bam}"
      input: fname
      output: join_path(f"data/bams/cond{condition}_{i}.bam")
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
      with open({output}, "w") as f:
        for region in {params.regions}:
          f.write(f"{region}\n")
