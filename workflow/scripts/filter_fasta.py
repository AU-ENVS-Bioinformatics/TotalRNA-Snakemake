from Bio import SeqIO

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

exclude_ids = snakemake.input.exclude
infasta = snakemake.input.fasta
outfile = snakemake.output[0]
# Read the exclude ids from file. There's one header per line.

exclude_ids = [line.strip() for line in open(exclude_ids, "r")]

acc = 0
with open(infasta, "r") as fasta, open(outfile, "w") as out:
    for record in SeqIO.parse(fasta, "fasta"):
        if record.id not in exclude_ids:
            SeqIO.write(record, out, "fasta")
        else:
            acc += 1
print(f"Removed {acc} sequences from {infasta}")
