from pathlib import Path
import sys
import pandas as pd

infile = snakemake.input[0]
outfile = snakemake.output[0]
evalue = snakemake.params.evalue

column_names = [
    "target_name",
    "accession",
    "query_name",
    "accession2",
    "mdl",
    "mdl_from",
    "mdl_to",
    "seq_from",
    "seq_to",
    "strand",
    "trunc",
    "pass",
    "gc",
    "bias",
    "score",
    "E-value",
    "inc",
    "description_length",
    "description_path",
]


with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f
    print(f"Reading file: {infile}", file=sys.stderr)
    df = pd.read_csv(
        infile,
        delimiter="\s+",
        header=None,
        comment="#",
        names=column_names,
    )
    # Join description_path column and description_length column into one
    df["description"] = df["description_path"].str.cat(
        df["description_length"].astype(str), sep=" "
    )
    # Drop description_path and description_length columns
    df = df.drop(columns=["description_path", "description_length"])
    print(f"Shape before: {df.shape}", file=sys.stderr)
    # Filter by E-value
    print(f"Filtering by E-value <= {evalue}", file=sys.stderr)
    df = df[df["E-value"] <= evalue]
    print(f"Shape before filtering by E-value: {df.shape}", file=sys.stderr)
    # Order by score and E-value
    df = df.sort_values(by=["score", "E-value"], ascending=[False, True])
    # Remove duplicates of column "target_name"
    df = df.drop_duplicates(subset=["target_name"])
    print(f"Shape after removing duplicates: {df.shape}", file=sys.stderr)
    # Write to file
    print(f"Writting to {outfile}", file=sys.stderr)
    df.to_csv(outfile, sep="\t", index=False)
