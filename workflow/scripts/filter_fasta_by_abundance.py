from pathlib import Path
import sys
import pandas as pd
from Bio import SeqIO

# Header column name
ID_COLUMN = "ContigID"
# Filepaths
filtered_table = pd.read_table(snakemake.input.table)       
infile = Path(str(snakemake.input.fasta))
outfile = Path(str(snakemake.output))
# Get desired headers
headers = filtered_table[ID_COLUMN].values

# Filter sequences
included = 0
with outfile.open('w') as out:
    for index, record in enumerate(SeqIO.parse(infile.open(), 'fasta')):
        if record.id in headers:
            SeqIO.write(record, out, "fasta")
            included+=1

# Print log
with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f
    print(f"{str(index)} sequences were processed")
    print(f"Only {str(included)} sequences were included")