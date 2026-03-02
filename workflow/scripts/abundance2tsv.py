import sys

import duckdb
import pandas as pd


database = snakemake.input.db
output_tsv = snakemake.output[0]
table_name = snakemake.params.table
log = snakemake.log

with open(log[0], "w") as f:
    sys.stderr = sys.stdout = f
    print(f"Exporting {table_name} from {database} to {output_tsv}")

    conn = duckdb.connect(database=database, read_only=True)
    try:
        df = conn.execute(
            f"SELECT contig AS ContigID, sample, mapped_reads FROM {table_name}"
        ).fetchdf()
    finally:
        conn.close()

    if df.empty:
        print("No rows found in table. Writing empty TSV with contig column.")
        pd.DataFrame(columns=["ContigID"]).to_csv(output_tsv, sep="\t", index=False)
        sys.exit(0)
    
    samples = sorted(df["sample"].unique().tolist())
    wide = df.pivot_table(
        index="ContigID",
        columns="sample",
        values="mapped_reads",
        aggfunc="sum",
        fill_value=0,
    )
    wide = wide.reindex(columns=samples)
    wide.reset_index(inplace=True)

    if samples:
        non_zero = wide[samples].sum(axis=1) > 0
        wide = wide.loc[non_zero]

    wide.to_csv(output_tsv, sep="\t", index=False)
    print(f"Done. Wrote {len(wide)} rows and {len(samples)} samples.")
