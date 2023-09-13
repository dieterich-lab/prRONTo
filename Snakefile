import pandas as pd
from snakemake.utils import validate
import yaml


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
  input: workflow.basedir + "/schemas/config.yaml"
  run:
    with open(input[0]) as f:
      schema = yaml.load(f, Loader=yaml.SafeLoader)
      print(yaml.dump(schema, default_flow_style=False))


rule describe_pep:
  input: workflow.basedir + "/schemas/pep.yaml"
  run:
    with open(input[0]) as f:
      schema = yaml.load(f, Loader=yaml.SafeLoader)
      print(yaml.dump(schema, default_flow_style=False))
