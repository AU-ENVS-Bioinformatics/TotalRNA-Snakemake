from pathlib import Path
from typing import List

from snakemake.io import expand


# =========================================================
# Kraken2 configuration helpers
# =========================================================

def kraken_database_names(config: dict) -> List[str]:
    """
    Return the names of all configured Kraken2 databases.

    Example:
        ["silva", "pr2"]
    """

    databases = config["rRNA"]["kraken2"].get("databases", {})

    if not databases:
        raise ValueError(
            "No Kraken2 databases were configured under "
            "rRNA -> kraken2 -> databases"
        )

    return list(databases.keys())


def get_kraken_database(database: str, config: dict) -> str:
    """
    Return the path belonging to a Kraken2 database name.
    """

    databases = config["rRNA"]["kraken2"]["databases"]

    if database not in databases:
        raise ValueError(
            f"Unknown Kraken2 database '{database}'. "
            f"Available databases: {', '.join(databases.keys())}"
        )

    database_path = databases[database].get("path")

    if not database_path:
        raise ValueError(
            f"No path was configured for Kraken2 database "
            f"'{database}'"
        )

    return database_path


def run_bracken(database: str, config: dict) -> bool:
    """
    Return whether Bracken should run for a database.
    """

    databases = config["rRNA"]["kraken2"]["databases"]

    if database not in databases:
        raise ValueError(
            f"Unknown Kraken2 database '{database}'"
        )

    return bool(
        databases[database].get("run_bracken", False)
    )


def bracken_database_names(config: dict) -> List[str]:
    """
    Return only databases for which Bracken is enabled.
    """

    return [
        database
        for database in kraken_database_names(config)
        if run_bracken(database, config)
    ]


# =========================================================
# Reports used for downstream processing
# =========================================================

def taxonomy_report(
    results_dir: Path,
    sample: str,
    database: str,
    config: dict,
) -> str:
    """
    Return the report used for downstream processing.

    When Bracken is enabled:
        sample.database.bracken.k2report

    When Bracken is disabled:
        sample.database.k2report
    """

    classification_dir = (
        Path(results_dir)
        / "rRNA"
        / sample
        / "classification"
        / database
    )

    if run_bracken(database, config):
        return str(
            classification_dir
            / f"{sample}.{database}.bracken.k2report"
        )

    return str(
        classification_dir
        / f"{sample}.{database}.k2report"
    )


def taxonomy_reports(
    results_dir: Path,
    samples: List[str],
    database: str,
    config: dict,
) -> List[str]:
    """
    Return all reports required to construct one combined
    BIOM table for one database.
    """

    return [
        taxonomy_report(
            results_dir=results_dir,
            sample=sample,
            database=database,
            config=config,
        )
        for sample in samples
    ]


# =========================================================
# Expected Kraken2 outputs
# =========================================================

def kraken_outputs(
    results_dir: Path,
    samples: List[str],
    config: dict,
) -> List[str]:
    """
    Return Kraken2 outputs for every sample and database.
    """

    databases = kraken_database_names(config)

    outputs: List[str] = []

    outputs += expand(
        (
            f"{results_dir}/rRNA/{{sample}}/classification/"
            "{database}/{sample}.{database}.k2report"
        ),
        sample=samples,
        database=databases,
    )

    outputs += expand(
        (
            f"{results_dir}/rRNA/{{sample}}/classification/"
            "{database}/{sample}.{database}.kraken"
        ),
        sample=samples,
        database=databases,
    )

    return outputs


# =========================================================
# Expected Bracken outputs
# =========================================================

def bracken_outputs(
    results_dir: Path,
    samples: List[str],
    config: dict,
) -> List[str]:
    """
    Return Bracken outputs only for databases where
    run_bracken is true.
    """

    databases = bracken_database_names(config)

    if not databases:
        return []

    outputs: List[str] = []

    outputs += expand(
        (
            f"{results_dir}/rRNA/{{sample}}/classification/"
            "{database}/{sample}.{database}.bracken.tsv"
        ),
        sample=samples,
        database=databases,
    )

    outputs += expand(
        (
            f"{results_dir}/rRNA/{{sample}}/classification/"
            "{database}/{sample}.{database}.bracken.k2report"
        ),
        sample=samples,
        database=databases,
    )

    return outputs


# =========================================================
# Expected downstream outputs
# =========================================================

def downstream_outputs(
    results_dir: Path,
    config: dict,
) -> List[str]:
    """
    Return downstream taxonomy outputs for every database.
    """

    databases = kraken_database_names(config)

    outputs: List[str] = []

    outputs += expand(
        (
            f"{results_dir}/rRNA/taxonomy/"
            "{database}/taxonomy.biom"
        ),
        database=databases,
    )

    outputs += expand(
        (
            f"{results_dir}/rRNA/taxonomy/"
            "{database}/taxonomy.tsv"
        ),
        database=databases,
    )

    outputs += expand(
        (
            f"{results_dir}/rRNA/taxonomy/"
            "{database}/taxonomy_qualitative.txt"
        ),
        database=databases,
    )

    outputs += expand(
        (
            f"{results_dir}/rRNA/taxonomy/"
            "{database}/taxonomy_observations.txt"
        ),
        database=databases,
    )

    outputs += expand(
        (
            f"{results_dir}/rRNA/taxonomy/"
            "{database}/phyloseq_raw.rds"
        ),
        database=databases,
    )

    outputs += expand(
        (
            f"{results_dir}/rRNA/taxonomy/"
            "{database}/phyloseq_filtered.rds"
        ),
        database=databases,
    )

    return outputs


# =========================================================
# Main rRNA target function
# =========================================================

def rrna_outputs(
    results_dir: Path,
    samples: List[str],
    config: dict,
) -> List[str]:
    """
    Return all required outputs for the rRNA module.
    """

    method = config["rRNA"]["method"]

    if method != "kraken2":
        raise ValueError(
            f"Unsupported rRNA classification method: {method}"
        )

    outputs: List[str] = []

    outputs += kraken_outputs(
        results_dir=results_dir,
        samples=samples,
        config=config,
    )

    outputs += bracken_outputs(
        results_dir=results_dir,
        samples=samples,
        config=config,
    )

    outputs += downstream_outputs(
        results_dir=results_dir,
        config=config,
    )

    return outputs