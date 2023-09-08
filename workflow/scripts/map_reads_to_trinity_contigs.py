from pathlib import Path
from snakemake.shell import shell

log = snakemake.log[0]

left_fasta_files = snakemake.input.left
right_fasta_files = snakemake.input.right
contigs = snakemake.input.contigs
threads = snakemake.threads
