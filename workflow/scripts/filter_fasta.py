from Bio import SeqIO
import pandas as pd

exclude_tsv = snakemake.input.exclude
infasta = snakemake.input.fasta
outfile = snakemake.output[0]

# Read target names from TSV
df = pd.read_csv(exclude_tsv, sep='\t')
exclude_ids = set(df['target_name'].tolist())

acc = 0
with open(infasta, "r") as fasta, open(outfile, "w") as out:
    for record in SeqIO.parse(fasta, "fasta"):
        if record.id not in exclude_ids:
            SeqIO.write(record, out, "fasta")
        else:
            acc += 1
print(f"Removed {acc} sequences from {infasta}")