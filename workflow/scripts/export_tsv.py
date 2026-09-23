import sys
import duckdb
import pandas as pd


database = snakemake.input.db
output_tsv = snakemake.output
log = snakemake.log

model = snakemake.params.model
threshold = snakemake.params.threshold

TABLES = ["mapped_reads", "read_length"]

    
def log_message(message):
    print(message, file=sys.stderr)


def export_sample_table(database_path, output_path, model="", threshold=0):
    table_name = "mapped_reads"
    conn = duckdb.connect(database=database_path)

    log_message(f"[{table_name}] Starting export from '{database_path}' to '{output_path}'")
    if model == "rpm":
        condition = "tr.total_mapped_reads / 1e6"
    elif model == "rpkm":
        condition = "tr.total_mapped_reads * cl.contig_length / 1e9"
    elif model == "rpkm_":
        condition = "tr.total_mapped_reads * cl.contig_length / rl.read_length  / 1e6"
    elif model == "tpm":
        condition = "cl.contig_length * td.denom / 1e6"
    else:
        condition = "1.0"

    query = f"""
    CREATE OR REPLACE VIEW {model}_{threshold} AS
    WITH total_reads AS (
    SELECT SUM(mapped_reads) AS total_mapped_reads
    FROM mapped_reads
    ),
    tpm_denominator AS (
        SELECT SUM(mr.mapped_reads * 1.0 / cl.contig_length) AS denom
        FROM mapped_reads mr
        JOIN contig_length cl
        ON mr.contig = cl.contig
    )
    SELECT
        mr.contig AS contig,
        mr.sample,
        mr.mapped_reads
    FROM mapped_reads mr
    JOIN read_length rl
        ON mr.contig = rl.contig
    JOIN contig_length cl
        ON mr.contig = cl.contig
    CROSS JOIN total_reads tr
    CROSS JOIN tpm_denominator td
    WHERE tr.total_mapped_reads > 0
    AND mr.mapped_reads >= {condition} * {threshold};
    """
    log_message(f"[{table_name}] Running query: {query}")

    conn.execute(query)
    df = conn.execute(f"SELECT * FROM {model}_{threshold}").fetchdf()
    log_message(f"[{table_name}] Retrieved {len(df)} rows")

    if df.empty:
        log_message(f"[{table_name}] No rows found. Writing empty TSV with 'contig' column")
        pd.DataFrame(columns=["contig"]).to_csv(output_path, sep="\t", index=False)
        log_message(f"[{table_name}] Finished with empty export")
        return

    log_message(f"[{table_name}] Pivoting table by sample")
    df = df.pivot(index="contig", columns="sample", values="mapped_reads").fillna(0)

    log_message(f"[{table_name}] Writing pivoted TSV")
    df.to_csv(output_path, sep="\t", index=True)
    log_message(
        f"[{table_name}] Done. Wrote {len(df)} contigs and {len(df.columns)} sample columns"
    )

    conn.close()

def export_contig_length(database_path, output_path, conn):
    table_name = "contig_length"
    log_message(f"[{table_name}] Starting export from '{database_path}' to '{output_path}'")
    query = "SELECT contig AS contig, contig_length FROM contig_length"
    log_message(f"[{table_name}] Running query: {query}")

    df = conn.execute(query).fetchdf()
    log_message(f"[{table_name}] Retrieved {len(df)} rows")

    if df.empty:
        log_message(f"[{table_name}] No rows found. Writing empty TSV with 'contig' column")
        pd.DataFrame(columns=["contig"]).to_csv(output_path, sep="\t", index=False)
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

    log_message("Opening DuckDB connection")

    try:
        export_sample_table(database, output_tsv[0], model=model, threshold=threshold)
        log_message("All exports completed successfully")
    finally:
        log_message("Closing DuckDB connection")
        log_message("TSV export script finished")
    

