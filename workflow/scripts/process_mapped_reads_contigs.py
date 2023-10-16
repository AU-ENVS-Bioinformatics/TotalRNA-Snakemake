from pathlib import Path
import sys
import pandas as pd

files = snakemake.input
log = snakemake.log[0]
params = snakemake.params

extension_fwd = params.extension_fwd
extension_rev = params.extension_rev


def regex_filename(filename):
    # Get filename
    filename = Path(filename).name
    # Remove everything after either extension_fwd or extension_rev
    filename = filename.split(extension_fwd)[0]
    filename = filename.split(extension_rev)[0]
    return filename


def read_samtools_idxstats(idxstats_file):
    df = pd.read_csv(idxstats_file, sep="\t", header=None, comment="*")
    df.columns = ["contig", "length", "mapped_reads", "unmapped_reads"]
    df["sample"] = regex_filename(idxstats_file)
    # return only sample, contig and mapped_reads
    return df[["sample", "contig", "mapped_reads"]]


with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f
    print("Reading samtools idxstats files...", file=sys.stderr)
    print(f"There are {len(files)} files", file=sys.stderr)
    # Read all samtools idxstats files into a single dataframe
    df = pd.concat([read_samtools_idxstats(f) for f in files], ignore_index=True)
    print("Done reading samtools idxstats files", file=sys.stderr)
    print(f"Table has shape: {df.shape}", file=sys.stderr)
    print("Grouping by sample and contig...", file=sys.stderr)
    # Group by sample and contig and sum mapped reads
    df = df.groupby(["sample", "contig"]).sum().reset_index()
    print("Done grouping by sample and contig", file=sys.stderr)
    print(f"Table has shape: {df.shape}", file=sys.stderr)
    # Pivot wider, with one contig per row the mapped_reads of each sample as a new column
    print("Pivoting wider...", file=sys.stderr)
    df = df.pivot(index="contig", columns="sample", values="mapped_reads").reset_index()
    # Rename contig to ContigID
    df = df.rename(columns={"contig": "ContigID"})
    print("Done pivoting wider", file=sys.stderr)
    print(f"Table has shape: {df.shape}", file=sys.stderr)
    # Write to file
    print("Writing to file...", file=sys.stderr)
    df.to_csv(snakemake.output[0], sep="\t", index=False)
    print("Done writing to file", file=sys.stderr)
