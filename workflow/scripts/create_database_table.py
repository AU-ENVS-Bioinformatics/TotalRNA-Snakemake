import sys
import duckdb
import re
import pandas as pd

database = snakemake.output.database
log = snakemake.log
single_sample = snakemake.input.single_sample
idxstats = snakemake.input.counts_length_tsv[0]
reads = snakemake.input.counts_length_tsv[-1]
counts_tsv = re.sub(r"GP\d{3}_(rev|fwd)", "*", idxstats)
length_tsv = re.sub(r"GP\d{3}_(rev|fwd)", "*", reads)
table_name = ["abundance", "total_read_length"]


value = ["mapped_reads", "read_length"]
input = [counts_tsv, length_tsv]

def write_to_table(database, table_name, value, input):
    query = f"""
    SELECT sample AS sample, contig AS contig, SUM({value}) as {value}
    FROM read_csv_auto('{input}')
    GROUP BY sample, contig
    HAVING SUM({value}) > 0"""
    print("Connecting to database and writing table...", file=sys.stderr)
    conn = duckdb.connect(database)
    conn.execute(f"INSERT INTO {table_name} {query}")
    conn.commit()
    return conn


def init_database_sample_contig(database, table_name, value):
    with duckdb.connect(database) as conn:
        conn.execute(f"""
            DROP TABLE IF EXISTS {table_name};
        """)
        conn.commit()
        conn.execute(f"""
        CREATE TABLE IF NOT EXISTS {table_name} (
            sample TEXT,
            contig TEXT,
            {value} INT
        )
    """)
        conn.commit()
        
def init_database_contig(database, table_name):
    with duckdb.connect(database) as conn:
        conn.execute(f"""
            DROP TABLE IF EXISTS {table_name};
        """)
        conn.commit()
        conn.execute(f"""
        CREATE TABLE IF NOT EXISTS {table_name} (
            contig TEXT,
            contig_length INT
        )
    """)
        conn.commit()

with open(log[0], "w") as f:
    sys.stderr = sys.stdout = f
    print("Processing mapped reads to contigs...", file=sys.stderr)

    for tbl, val, inp in zip(table_name, value, input):
        init_database_sample_contig(database, tbl, val)
        conn = write_to_table(database, tbl, val, inp)
        print(f"Done writing table {tbl}", file=sys.stderr)
        conn.close()
    print("Finished processing mapped reads to contigs.", file=sys.stderr)
    print("Processing contig lengths...", file=sys.stderr)
    init_database_contig(database, "contig_length")
    df = pd.read_csv(single_sample, sep="\t", header=None, usecols=[0, 1], names=["contig", "contig_length"])
    df.to_sql("contig_length", duckdb.connect(database), if_exists="replace", index=False)
    print("Finished writing contig_length table.", file=sys.stderr)