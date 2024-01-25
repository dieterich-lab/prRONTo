# options activate conda evn

# run pipeline
snakemake -c 4 --snakefile ../../Snakefile --config pepfile=data.yml --configfile=analysis.yml

