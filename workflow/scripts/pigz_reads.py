import sys
import os

infiles = snakemake.input.R1 + snakemake.input.R2
outfiles = snakemake.output.R1 + snakemake.output.R2
threads = snakemake.threads

with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f
    for old, new in zip(infiles, outfiles):
        cmd = f"pigz -dkf -p{threads} < {old} > {new}"
        os.system(cmd)
        print(cmd)
        