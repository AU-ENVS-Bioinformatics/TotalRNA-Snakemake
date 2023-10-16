import pandas as pd
import sys


def parse_cmsearch_line(line):
    line = line.split()
    return pd.Series(
        {
            "target_name": line[0],
            "query_name": line[2],
            "accession": line[3],
            "score": float(line[14]),
            "evalue": float(line[15]),
        }
    )


infile = snakemake.input[0]
outfile = snakemake.output[0]
evalue = snakemake.params.evalue

with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f
    print(f"Reading file: {infile}", file=sys.stderr)
    ## Create dataframe from a list of dictionaries

    df = pd.DataFrame(
        [
            parse_cmsearch_line(line)
            for line in open(infile)
            if not line.startswith("#")
        ],
    )
    print(f"Shape before: {df.shape}", file=sys.stderr)
    # Filter by E-value
    print(f"Filtering by E-value <= {evalue}", file=sys.stderr)
    df = df[df["evalue"] <= evalue]
    print(f"Shape before filtering by E-value: {df.shape}", file=sys.stderr)
    # Order by score and E-value
    df = df.sort_values(by=["score", "evalue"], ascending=[False, True])
    # Remove duplicates of column "target_name"
    df = df.drop_duplicates(subset=["target_name"])
    print(f"Shape after removing duplicates: {df.shape}", file=sys.stderr)
    # Write to file
    print(f"Writting to {outfile}", file=sys.stderr)
    df.to_csv(outfile, sep="\t", index=False)
