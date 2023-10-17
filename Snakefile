import pandas as pd
from snakemake.utils import validate
import yaml


validate(config, "schemas/config.yaml")
pepfile: config["pepfile"]
pepschema: "schemas/pep.yaml"


shell.prefix(f"export PRONTO_DIR={workflow.basedir} ; ")


include: "rules/common.smk"
include: "rules/samtools.smk"
include: "rules/jacusa2.smk"
include: "rules/meta.smk"
include: "rules/plot.smk"
include: "rules/report.smk"


rule all:
  input: auto_targets()


rule original:
  input: original_targets()


rule downsampling:
  input: downsampling_targets()


def describe_yaml(fname):
    with open(input[0]) as f:
      schema = yaml.load(f, Loader=yaml.SafeLoader)
      print(yaml.dump(schema, default_flow_style=False))


rule describe_config:
  input: workflow.basedir + "/schemas/config.yaml"
  run:
    describe_yaml(input[0])


rule describe_pep:
  input: workflow.basedir + "/schemas/pep.yaml"
  run:
    describe_yaml(input[0])
