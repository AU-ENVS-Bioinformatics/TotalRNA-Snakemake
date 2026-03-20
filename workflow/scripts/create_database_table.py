import sys
import duckdb

database = snakemake.output.database
log = snakemake.log
single_sample = snakemake.input.idxstats[0]
dir_path = snakemake.params.dir_path

INPUT_PATTERNS = f"{dir_path}/*_sorted.bam.idxstats"


def log_message(message):
    print(message, file=sys.stderr)


def table_read_csv(input_pattern):
    return f"""
        FROM read_csv(
        '{input_pattern}',
        delim='\t',
        header=False,
        columns={{
            'contig': 'TEXT',
            'contig_length': 'INTEGER',
            'mapped': 'INTEGER',
            'unmapped': 'INTEGER'
        }},
        ignore_errors=True,
        filename=True
    )
    """


def write_contig_length_table(conn, bam_file):
    log_message("[contig_length] Building insert query")
    query = f"""
    SELECT contig, contig_length
    {table_read_csv(bam_file)}
    WHERE contig_length > 0
    GROUP BY contig, contig_length
    """
    log_message(f"[contig_length] Inserting rows from pattern: {bam_file}")
    conn.execute(f"INSERT INTO contig_length {query}")
    log_message("[contig_length] Insert completed")


def write_mapped_reads_table(conn, mapped_pattern):
    log_message("[mapped_reads] Building insert query")
    query = f"""
    SELECT contig,regexp_replace(
    regexp_extract(filename, '[^/]+$'),
    '_sorted\\.bam\\.idxstats$', '') AS sample, mapped AS mapped_reads
    {table_read_csv(mapped_pattern)}
    WHERE mapped > 0
    GROUP BY contig, sample, mapped_reads
    """
    log_message(f"[mapped_reads] Inserting rows from pattern: {mapped_pattern}")
    conn.execute(f"INSERT INTO mapped_reads {query}")
    log_message("[mapped_reads] Insert completed")


def init_database(database):
    log_message(f"Initializing DuckDB database: {database}")
    conn = duckdb.connect(database)

    log_message("Dropping old tables if they exist")
    conn.execute("DROP TABLE IF EXISTS mapped_reads")
    conn.execute("DROP TABLE IF EXISTS read_length")
    conn.execute("DROP TABLE IF EXISTS contig_length")

    log_message("Creating table: contig_length")
    conn.execute("""
    CREATE TABLE contig_length (
        contig TEXT,
        contig_length INT
    )
    """)

    log_message("Creating table: mapped_reads")
    conn.execute("""
    CREATE TABLE mapped_reads (
        contig TEXT,
        sample TEXT,
        mapped_reads INT
    )
    """)

    log_message("Creating table: read_length")
    conn.execute("""
    CREATE TABLE read_length (
        contig TEXT,
        sample TEXT,
        read_length INT
    )
    """)

    log_message("Database initialization complete")

    return conn


with open(log[0], "w") as f:
    sys.stderr = sys.stdout = f
    log_message("Starting database table creation script")
    log_message(f"Output database: {database}")
    log_message(f"Input contig pattern: {single_sample}")
    log_message(f"Input mapped pattern: {INPUT_PATTERNS}")
    log_message(f"Log file: {log[0]}")

    conn = None
    try:
        log_message("Opening database connection and preparing schema")
        conn = init_database(database)

        log_message("Writing contig_length data")
        write_contig_length_table(conn, single_sample)

        log_message("Writing mapped_reads data")
        write_mapped_reads_table(conn, INPUT_PATTERNS)

        log_message("Committing transaction")
        conn.commit()
        log_message("Commit completed successfully")
    finally:
        if conn is not None:
            log_message("Closing database connection")
            conn.close()
        log_message("Database table creation script finished")