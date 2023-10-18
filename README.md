prRONTo - precise rRNA modification analysis on ONTechnoligies
==============================================================

[prRONTo](https://github.com/dieterich-lab/prRONTo) is a snakemake workflow 
for comprehensive identification of diverse ribosomal RNA modifications
 by targeted nanopore direct RNA sequencing and JACUSA2.

See original publication for details [DOI: 10.1080/15476286.2023.2248752](https://www.tandfonline.com/doi/full/10.1080/15476286.2023.2248752)

# Installation

Chekout the repository and install required packages (see `conda.yaml`).

```
git clone https://github.com/dieterich-lab/prRONTo
cd prRONTo
```

We recommend to use [Conda](https://conda.io) to install the required packages.

## Requirements

Internally, [prRONTo](https://github.com/dieterich-lab/prRONTo) uses 
[JACUSA2](https://dieterich-lab/JACUSA2) and 
[JACUSA2helper](https://github.com/dieterich-lab/JACUSA2helper) to
detect RNA modifications and a collection of 
[R](https://www.r-project.org) and [Python](https://www.python.org) scripts for processing.

Required dependencies can be installed via [Conda](https://conda.io):
```
conda env create -n pronto -f conda.yaml
conda activate pronto
```

# Usage

A a sample description via a [PEP](https://pep.databio.org/en/2.0.0) file and a sample table is required.
Furthermore, the reference FASTA and a customized RNA modification file corresponding to the sequencing data is required.

```
snakemake -c 1 -f <SNAKEFILE> --config pep=<PEPFILE> [--configfile=<CONFIG_FILE>]
```
The `<CONFIG_FILE>` defines and adjusts parameters of the analysis whereas the `<PEPFILE>` is entirely sample specific.

## PEP file

The [PEP](https://pep.databio.org/en/2.0.0/) is in [yaml](https://yaml.org/) format.
In the following, a descriptive example is presented:

```yaml
pep_version: 2.0.0
sample_table: <SAMPLE_TABLE>      # REQUIRED, path to sample table

project: "Example"                # OPTIONAL, plain text, no special characters

pronto:
  regions: ["region1", "region2"] # REQUIRED, list of seqnames to scan for modifications
  output: <OUTPUT_DIRECTOR>       # REQUIRED, path where results will be written to
  mods: <MODS_FILE>               # REQUIRED, path to existing modification annotation
  ref: <REF_FASTA>                # REQUIRED, path to reference FASTA

  # values for condition 1 and 2 
  # must exist in
  # the column "condition" in the sample table!
  condition1: "condition1"        # REQUIRED, value from col. "condition" in sample table
  condition2: "condition2"        # REQUIRED, value from col. "condition", in sample table
```

See `example/human/pep.yaml` for an example.

## Sample table

A minimal sample table is provided in the following:

| sample_name | condition  | bam       |
| ------------| ---------- | --------- |
| sample_1_1  | condition1 | <BAM_1_1> |
| sample_1_2  | condition1 | <BAM_1_2> |
| sample_2_1  | condition2 | <BAM_2_1> |
| sample_2_2  | condition2 | <BAM_2_2> |

See `example/human/sample_table.yaml` for an example.

## Analysis config

In the following, a descriptive example is presented:

```yaml
downsampling:   # OPTIONAL
  # target coverage
  reads: [1000, ]
  # seeds for read sampling with samtools
  # the number of seeds corresponds to the number of repetitions
  seed: ["CfdCY", "gJm9e", "CZ8X7", "1Cdq4", "6pod1", "RtDnQ", "AdOWe", ]
jacusa2:
  features: ["M", "MDI"]    # Features to use for the LOF analysis
# LOF specific parameters
lof:
  - neighbors: 20
    contamination: 0.001
  - neighbors: 20
    contamination: 0.002
```

# Modification file format

A tab-delimeted file with three columns: `seq_id`, `pos`(0-indexed), and `mod`.

In the following a section from `examples/data/human_rrna.tsv`:

| seq_id             | pos | mod |
| ------------------ | --- | --- |
| NR_003286_RNA18SN5 | 26  | Am  |
| NR_003286_RNA18SN5 | 33  | psU |
| NR_003286_RNA18SN5 | 35  | psU |
| NR_003286_RNA18SN5 | 92  | psU |
| ...                | ... | ... |

# Report

The report is organised in to 4 sections:

* Results
* Read
* Modifications and 
* Session.

## Results

The results section contains the putative RNA modifications as outliers stratified 
by region and utilized feature.
Additionaly, known modifications, optional downsampling info and LOF parameters are presented.

## Reads

This section is mainly based on `samtools stats|coverage`. 
Properties of reads and summary statistics are presented:

* total reads,
* mapped reads,
* read length,
* mapping quaylity, 
* ...

## Modifications

Summary of provided RNA modification file.

## Session

Paths to config files and installed software can be checked in this section.

# Example

In the following, a chain of commands to run a toy example:

1. Clone repository: `git clone https://github.com/dieterich-lab/prRONTo`
2. `cd prRONTo`
3. Install conda env.: `conda env create -n pronto -f conda.yaml`
4. Activate conda env.: `conda activate pronto`
5. Run prRONTo: `cd example/human ; run_exampe.sh``
6. Open `output/report/report.html`with your favorite browser and check the results.
