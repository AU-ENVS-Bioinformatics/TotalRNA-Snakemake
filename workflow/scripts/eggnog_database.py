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

KO_file = output.KO_file
gene_file = output.gene_file
eggnog_output_file = output.eggnog_output_file

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
            skiprows=4,
            header=0,
            skipfooter=3,
        )

        eggnog.rename(
            columns={
                "#query": "gene",
                "max_annot_lvl": "taxonomy",
                "Description": "function",
            },
            inplace=True,
        )

        # Normalize gene names
        eggnog["gene"] = eggnog["gene"].map(normalize_gene)

        # Filter by quality thresholds
        eggnog = eggnog.query("evalue <= 1e-10 and score >= 80")

        return eggnog
    except Exception as e:
        print(f"Error loading eggNOG file: {e}", file=sys.stderr)
        raise


def main():
    """Main execution function."""
    with open(log[0], "w") as f:
        sys.stderr = sys.stdout = f
        try:
            # Connect to database
            con = duckdb.connect(DB)
            print(f"Connected to database: {DB}", file=sys.stderr)

            # Load and filter eggNOG annotations
            eggnog_df = load_and_filter_eggnog(EGGNOG)
            print(
                f"Loaded {len(eggnog_df)} filtered eggNOG annotations", file=sys.stderr
            )

            # Populate eggnog_output table
            con.execute(
                "CREATE OR REPLACE TABLE eggnog_output AS SELECT * FROM eggnog_df"
            )
            print("eggNOG annotations loaded to database", file=sys.stderr)

            con.execute(
                """DROP VIEW IF EXISTS eggnog_best; CREATE VIEW eggnog_best AS
                        SELECT gene, function, Preferred_name, KEGG_ko, taxonomy
            FROM (
                SELECT *,
                       ROW_NUMBER() OVER (
                           PARTITION BY gene
                           ORDER BY evalue ASC, score DESC
                       ) rn
                FROM eggnog_output
            )
            WHERE rn = 1"""
            )

            ko_gene_data = con.execute(
                f"""
                SELECT
                    c.gene,
                    c.sample,
                    c.mapped_reads,
                    COALESCE(e.KEGG_ko, 'NA') AS KEGG_ko,
                    COALESCE(e.function, 'NA') AS function,
                    COALESCE(e.Preferred_name, 'NA') AS Preferred_name,
                    COALESCE(e.taxonomy, 'NA') AS taxonomy
                FROM (
                    SELECT
                        REPLACE(contig, substr(contig, instr(contig, '_i')), '') AS gene,
                        sample,
                        mapped_reads
                    FROM mapped_reads
                ) c
                LEFT JOIN eggnog_best e
                ON c.gene = e.gene
                WHERE e.KEGG_ko IS NOT NULL AND e.KEGG_ko NOT IN ('NA', '-', '')
                """
            ).fetchdf()

            df = ko_gene_data.pivot_table(
                index="KEGG_ko",
                columns="sample",
                values="mapped_reads",
                aggfunc="sum",
                fill_value=0,
            ).reset_index()
            df.to_csv(KO_file, sep="\t", index=False)

            df = ko_gene_data.pivot_table(
                index="gene",
                columns="sample",
                values="mapped_reads",
                aggfunc="sum",
                fill_value=0,
            ).reset_index()
            df.to_csv(gene_file, sep="\t", index=False)
            print("Gene and KO counts aggregated successfully", file=sys.stderr)

            df = (
                ko_gene_data.groupby(
                    ["gene", "function", "Preferred_name", "KEGG_ko", "taxonomy"]
                )
                .first()
                .reset_index()
                .drop(columns=["sample", "mapped_reads"])
            )
            df.to_csv(eggnog_output_file, sep="\t", index=False)

            con.close()

        except Exception as e:
            print(f"Error processing eggNOG data: {e}", file=sys.stderr)
            sys.exit(1)


if __name__ == "__main__":
    main()
