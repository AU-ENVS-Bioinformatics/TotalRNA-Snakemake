"""
Process eggNOG annotations and populate DuckDB database with gene and KO counts.
"""
import pandas as pd
import duckdb
import re
import sys
from typing import Optional

# Snakemake inputs and outputs
input = snakemake.input
log = snakemake.log
params = snakemake.params
output = snakemake.output

EGGNOG = input.eggnog
DB = input.database


def normalize_gene(gene_name: str) -> str:
    """
    Normalize gene names by removing isoform and protein suffixes.
    
    Args:
        gene_name: Raw gene identifier
        
    Returns:
        Normalized gene identifier
    """
    gene_name = re.sub(r"\.p\d+$", "", gene_name)
    gene_name = re.sub(r"_i\d+$", "", gene_name)
    return gene_name


def load_and_filter_eggnog(eggnog_file: str) -> pd.DataFrame:
    """
    Load eggNOG annotation file and filter by quality thresholds.
    
    Args:
        eggnog_file: Path to eggNOG annotation file
        
    Returns:
        Filtered DataFrame with eggNOG annotations
    """
    try:
        eggnog = pd.read_csv(
            eggnog_file,
            sep="\t",
            comment="#",
            header=None,
            usecols=[0, 2, 3, 5, 6, 7, 8, 10, 11],
            names=["gene", "evalue", "score", "taxonomy", "cog_category", 
                   "function", "preferred_name", "ec", "kegg_ko"]
        )
        
        # Normalize gene names
        eggnog["gene"] = eggnog["gene"].map(normalize_gene)
        
        # Filter by quality thresholds
        eggnog = eggnog.query("evalue <= 1e-10 and score >= 80")
        
        return eggnog
    except Exception as e:
        print(f"Error loading eggNOG file: {e}", file=sys.stderr)
        raise


def initialize_database(con: duckdb.DuckDBPyConnection) -> None:
    """
    Create database tables and clear existing data.
    
    Args:
        con: DuckDB database connection
    """
    con.execute(
        """
        CREATE TABLE IF NOT EXISTS deseq2_KO_counts (
            kegg_ko VARCHAR,
            sample VARCHAR,
            mapped_reads INTEGER
        );
        """
    )
    con.execute(
        """
        CREATE TABLE IF NOT EXISTS deseq2_gene_counts (
            gene VARCHAR,
            sample VARCHAR,
            mapped_reads INTEGER
        );
        """
    )
    con.execute("DELETE FROM deseq2_KO_counts")
    con.execute("DELETE FROM deseq2_gene_counts")


def populate_eggnog_table(eggnog_df: pd.DataFrame, con: duckdb.DuckDBPyConnection) -> None:
    """
    Populate the eggnog_output table with filtered annotations.
    
    Args:
        eggnog_df: Filtered eggNOG DataFrame
        con: DuckDB database connection
    """
    con.register("eggnog_input", eggnog_df)
    con.execute("CREATE OR REPLACE TABLE eggnog_output AS SELECT * FROM eggnog_input")

def aggregate_counts(con: duckdb.DuckDBPyConnection) -> None:
    """
    Aggregate gene and KO counts from abundance data and eggNOG annotations.
    
    Creates a temporary table joining abundance data with best eggNOG matches,
    then populates the final count tables.
    
    Args:
        con: DuckDB database connection
    """
    con.execute("""
        DROP TABLE IF EXISTS final_long;
    """)
    
    con.execute("""
        CREATE TEMP TABLE final_long AS
        WITH
        counts_gene AS (
            SELECT
                REPLACE(a.contig, substr(a.contig, instr(a.contig, '_i')), '') AS gene,
                a.sample,
                SUM(a.mapped_reads) AS mapped_reads
            FROM abundance a
            GROUP BY gene, sample
        ),
        eggnog_best AS (
            SELECT gene, kegg_ko
            FROM (
                SELECT *,
                       ROW_NUMBER() OVER (
                           PARTITION BY gene
                           ORDER BY evalue ASC, score DESC
                       ) rn
                FROM eggnog_output
            )
            WHERE rn = 1
        )
        SELECT
            c.gene,
            c.sample,
            c.mapped_reads,
            COALESCE(e.kegg_ko, 'NA') AS kegg_ko
        FROM counts_gene c
        LEFT JOIN eggnog_best e
        ON c.gene = e.gene;
    """)

    con.execute("""
        INSERT INTO deseq2_gene_counts (gene, sample, mapped_reads)
        SELECT gene, sample, mapped_reads
        FROM final_long;
    """)

    con.execute("""
        INSERT INTO deseq2_KO_counts (kegg_ko, sample, mapped_reads)
        SELECT kegg_ko, sample, SUM(mapped_reads) AS mapped_reads
        FROM final_long
        WHERE kegg_ko NOT IN ('NA', '-', '')
        GROUP BY kegg_ko, sample;
    """)


def main():
    """Main execution function."""
    with open(log[0], "w") as f:
        sys.stderr = sys.stdout = f
        try:
            # Connect to database
            con = duckdb.connect(DB)
            print(f"Connected to database: {DB}", file=sys.stderr)
            
            # Initialize database schema
            initialize_database(con)
            print("Database tables initialized", file=sys.stderr)
            
            # Load and filter eggNOG annotations
            eggnog_df = load_and_filter_eggnog(EGGNOG)
            print(f"Loaded {len(eggnog_df)} filtered eggNOG annotations", file=sys.stderr)
            
            # Populate eggnog_output table
            populate_eggnog_table(eggnog_df, con)
            print("eggNOG annotations loaded to database", file=sys.stderr)
            
            # Aggregate counts and populate final tables
            aggregate_counts(con)
            print("Gene and KO counts aggregated successfully", file=sys.stderr)
            
            con.close()
            
        except Exception as e:
            print(f"Error processing eggNOG data: {e}", file=sys.stderr)
            sys.exit(1)


if __name__ == "__main__":
    main()