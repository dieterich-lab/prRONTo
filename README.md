prRONTo - precise rRNA modification analysis on ONTechnoligies
=============================================================

# Installation
```
git clone https://github.com/dieterich-lab/prRONTo
cd prRONTo
```
## Requirements
Rquired dependencies can be installed via conda:
```
conda env create -n pronto -f conda.yaml
conda activate pronto
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

# Usage
```
snakemake -c 1 -f <SNAKEFILE> --pepfile <PEPFILE> (analysis|downsampling)
```

## Options
Use the following commands to get a description of available options for `config` and `pep`:
```
snakemake -c 1 -f <SNAKEFILE> describe_pep
snakemake -c 1 -f <SNAKEFILE> describe_config
```
