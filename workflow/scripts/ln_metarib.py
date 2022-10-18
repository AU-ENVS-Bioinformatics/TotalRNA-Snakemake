import os
infiles = snakemake.input.R1 + snakemake.input.R2
outfiles = snakemake.output.R1 + snakemake.output.R2
for old, new in zip(infiles, outfiles):
    os.system(f"cp {old} {new}.gz")
    os.system(f"gzip -d {new}.gz")
