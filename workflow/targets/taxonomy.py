from pathlib import Path
from typing import List

from snakemake.io import expand


# =========================================================
# SSU Kraken2 database helpers
# =========================================================

def kraken_database_names(
    config: dict,
) -> List[str]:
    """
    Return all configured SSU Kraken2 database names.

    Example:
        ["silva", "pr2"]
    """

    databases = config["taxonomy"]["SSU"]["kraken2"].get(
        "databases",
        {},
    )

    if not databases:
        raise ValueError(
            "No SSU Kraken2 databases configured under "
            "taxonomy -> SSU -> kraken2 -> databases"
        )

    return list(databases.keys())


def get_kraken_database(
    database: str,
    config: dict,
) -> str:
    """
    Return the path for one SSU Kraken2 database.
    """

    databases = config["taxonomy"]["SSU"]["kraken2"][
        "databases"
    ]

    if database not in databases:
        raise ValueError(
            f"Unknown SSU Kraken2 database '{database}'. "
            "Available databases: "
            + ", ".join(databases.keys())
        )

    database_path = databases[database].get("path")

    if not database_path:
        raise ValueError(
            f"No path configured for SSU Kraken2 database "
            f"'{database}'"
        )

    return database_path


# =========================================================
# Bracken reports used by Kraken-BIOM
# =========================================================

def taxonomy_report(
    results_dir: Path,
    sample: str,
    database: str,
) -> str:
    """
    Return the Bracken-adjusted SSU report used by Kraken-BIOM.
    """

    return str(
        Path(results_dir)
        / "Taxonomy_Profiling"
        / sample
        / "SSU"
        / database
        / "bracken"
        / f"{sample}.{database}.bracken.k2report"
    )


def taxonomy_reports(
    results_dir: Path,
    samples: List[str],
    database: str,
) -> List[str]:
    """
    Return all Bracken-adjusted SSU reports for one database.
    """

    return [
        taxonomy_report(
            results_dir=results_dir,
            sample=sample,
            database=database,
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
    Return SSU Kraken2 outputs for every sample and database.
    """

    databases = kraken_database_names(config)

    outputs: List[str] = []

    outputs += expand(
        f"{results_dir}/Taxonomy_Profiling/{{sample}}/SSU/{{database}}/kraken/{{sample}}.{{database}}.k2report",
        sample=samples,
        database=databases,
    )

    outputs += expand(
        f"{results_dir}/Taxonomy_Profiling/{{sample}}/SSU/{{database}}/kraken/{{sample}}.{{database}}.kraken",
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
    Return SSU Bracken outputs for every sample and database.
    """

    databases = kraken_database_names(config)

    outputs: List[str] = []

    outputs += expand(
        f"{results_dir}/Taxonomy_Profiling/{{sample}}/SSU/{{database}}/bracken/{{sample}}.{{database}}.bracken.tsv",
        sample=samples,
        database=databases,
    )

    outputs += expand(
        f"{results_dir}/Taxonomy_Profiling/{{sample}}/SSU/{{database}}/bracken/{{sample}}.{{database}}.bracken.k2report",
        sample=samples,
        database=databases,
    )

    return outputs


# =========================================================
# Expected SSU summary outputs
# =========================================================

def downstream_outputs(
    results_dir: Path,
    config: dict,
) -> List[str]:
    """
    Return combined SSU taxonomy outputs for every database.
    """

    databases = kraken_database_names(config)

    outputs: List[str] = []

    # taxonomy files created from kraken2/bracken reports
    outputs += expand(
        f"{results_dir}/Taxonomy_Profiling/summary/SSU/{{database}}/taxonomy/{{filename}}",
        database=databases,
        filename=[
            "taxonomy.biom",
            "taxonomy.tsv",
        ],
    )

    # statistics files created from biom summarize-table based on the kraken/bracken taxonomy files
    outputs += expand(
        f"{results_dir}/Taxonomy_Profiling/summary/SSU/{{database}}/statistics/{{filename}}",
        database=databases,
        filename=[
            "taxonomy_qualitative.txt",
            "taxonomy_observations.txt",
        ],
    )

    # converted taxonomy to phyloseq objects
    outputs += expand(
        f"{results_dir}/Taxonomy_Profiling/summary/SSU/{{database}}/phyloseq/{{filename}}",
        database=databases,
        filename=[
            "phyloseq_raw.rds",
            "phyloseq_filtered.rds",
        ],
    )

    return outputs


# =========================================================
# Main taxonomy output function
# =========================================================

def taxonomy_outputs(
    results_dir: Path,
    samples: List[str],
    config: dict,
) -> List[str]:
    """
    Return all required SSU taxonomy outputs.

    Workflow:
        SSU → Kraken2 → Bracken → BIOM → phyloseq
    """

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