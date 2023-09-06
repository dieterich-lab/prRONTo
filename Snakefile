import pandas as pd
from snakemake.utils import validate

validate(config, "schemas/config.yaml")
pepfile: config["pepfile"]
pepschema: "schemas/pep.yaml"


include: "rules/common.smk"
include: "rules/samtools.smk"
include: "rules/jacusa2.smk"
include: "rules/meta.smk"


rule all:
  input: auto_targets()


rule analysis:
  input: analysis_targets()


rule downsampling:
  input: downsampling_targets()


rule mixing:
  input: mixing_targets()


rule describe_config:
  input: "{workflow.basedir}/schemas/config.yaml"
  run:
    pass # TODO


rule describe_pep:
  input: "{workflow.basedir}/schemas/pop.yaml"
  run:
    pass # TODO
