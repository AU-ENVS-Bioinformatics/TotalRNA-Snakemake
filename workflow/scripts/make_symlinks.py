import sys
import os
from pathlib import Path

infiles = snakemake.input.R1 + snakemake.input.R2
outfiles = snakemake.output.R1 + snakemake.output.R2

with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f
    for src, dest in zip(infiles, outfiles):
        src, dest = Path(src).absolute(), Path(dest).absolute()
        dest.parent.mkdir(parents=True, exist_ok=True)
        print(src)
        print(dest)
        dest.symlink_to(src)
