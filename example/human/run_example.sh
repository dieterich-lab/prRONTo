# options activate conda evn

# run pipeline
snakemake -c 2 --snakefile ../../Snakefile --config pepfile=data.yml --configfile=analysis.yml

