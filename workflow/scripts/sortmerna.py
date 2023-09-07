from pathlib import Path
from snakemake.shell import shell

log = snakemake.log[0]

input_fasta = snakemake.input.fasta
database = snakemake.input.database
database_index = snakemake.input.database_index
aligned_prefix = snakemake.output.aligned[0][:-10]
not_aligned_prefix = snakemake.output.not_aligned[0][:-10]
workdir_prefix = Path(aligned_prefix).parent.parent / Path(aligned_prefix).name

shell(
    "sortmerna -ref {database} "
    "--idx-dir {database_index} "
    "--threads {snakemake.threads} "
    "--workdir {workdir_prefix} "
    "{snakemake.params.extra} --log "
    "--aligned {aligned_prefix} "
    "--other {not_aligned_prefix} "
    "--reads {input_fasta[0]} --reads {input_fasta[1]} "
    "2> {snakemake.output.stats} 1>&2"
)