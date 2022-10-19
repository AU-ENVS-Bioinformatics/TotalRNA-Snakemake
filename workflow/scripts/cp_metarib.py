import os
infiles = snakemake.input.R1 + snakemake.input.R2
outfiles = snakemake.output.R1 + snakemake.output.R2
threads = snakemake.threads
for old, new in zip(infiles, outfiles):
    os.system(f"pigz -dk -p{threads} {old} > {new}")