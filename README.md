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

# Usage
```
snakemake -c 1 -f <SNAKEFILE> --pepfile <PEPFILE> (analysis|downsampling|mixing)
```

## Options
Use the following commands to get a description of available options for `config` and `pep`:
```
snakemake -c 1 -f <SNAKEFILE> describe_pep
snakemake -c 1 -f <SNAKEFILE> describe_config
```
