import sys
import duckdb
import pandas as pd


database = snakemake.input.db
output_tsv = snakemake.output
log = snakemake.log

TABLES = ["mapped_reads", "read_length"]


def log_message(message):
    print(message, file=sys.stderr)


def export_sample_table(database_path, output_path, table_name, conn):
    log_message(f"[{table_name}] Starting export from '{database_path}' to '{output_path}'")
    query = f"SELECT contig AS ContigID, sample, {table_name} FROM {table_name}"
    log_message(f"[{table_name}] Running query: {query}")

    df = conn.execute(query).fetchdf()
    log_message(f"[{table_name}] Retrieved {len(df)} rows")

    if df.empty:
        log_message(f"[{table_name}] No rows found. Writing empty TSV with 'ContigID' column")
        pd.DataFrame(columns=["ContigID"]).to_csv(output_path, sep="\t", index=False)
        log_message(f"[{table_name}] Finished with empty export")
        return

    log_message(f"[{table_name}] Pivoting table by sample")
    df = df.pivot(index="ContigID", columns="sample", values=table_name).fillna(0)

    log_message(f"[{table_name}] Writing pivoted TSV")
    df.to_csv(output_path, sep="\t", index=True)
    log_message(
        f"[{table_name}] Done. Wrote {len(df)} contigs and {len(df.columns)} sample columns"
    )


def export_contig_length(database_path, output_path, conn):
    table_name = "contig_length"
    log_message(f"[{table_name}] Starting export from '{database_path}' to '{output_path}'")
    query = "SELECT contig AS ContigID, contig_length FROM contig_length"
    log_message(f"[{table_name}] Running query: {query}")

    df = conn.execute(query).fetchdf()
    log_message(f"[{table_name}] Retrieved {len(df)} rows")

    if df.empty:
        log_message(f"[{table_name}] No rows found. Writing empty TSV with 'ContigID' column")
        pd.DataFrame(columns=["ContigID"]).to_csv(output_path, sep="\t", index=False)
        log_message(f"[{table_name}] Finished with empty export")
        return

    log_message(f"[{table_name}] Writing TSV")
    df.to_csv(output_path, sep="\t", index=False)
    log_message(f"[{table_name}] Done. Wrote {len(df)} rows")

with open(log[0], "w") as f:
    sys.stderr = sys.stdout = f
    log_message("Starting TSV export script")
    log_message(f"Input database: {database}")
    log_message(f"Output files: {list(output_tsv)}")
    log_message(f"Log file: {log[0]}")

    log_message("Opening read-only DuckDB connection")
    conn = duckdb.connect(database=database, read_only=True)

    try:
        export_sample_table(database, output_tsv[0], TABLES[0], conn)
        export_sample_table(database, output_tsv[1], TABLES[1], conn)
        export_contig_length(database, output_tsv[2], conn)
        log_message("All exports completed successfully")
    finally:
        log_message("Closing DuckDB connection")
        conn.close()
        log_message("TSV export script finished")
    

