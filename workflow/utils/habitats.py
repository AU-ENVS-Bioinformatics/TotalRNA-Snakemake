"""
Helpers for habitat-specific taxonomy.

Provides habitat selection, configuration validation, database lookup,
Kraken2/Bracken settings, report paths, and paired-read selection.

Expected completion outputs remain in targets/habitats.py.
"""

from pathlib import Path
from typing import Any, List, Optional


# --------------------------
# Habitat selection
# --------------------------

VALID_HABITATS = {"Soil", "Glacier", "none"}


def selected_habitat(config: dict) -> Optional[str]:
    """Return the selected habitat, or None when disabled."""
    habitat = config.get("Habitat_refinement", {}).get(
        "habitat", "none"
    )

    if habitat not in VALID_HABITATS:
        raise ValueError(
            f"Unsupported Habitat_refinement.habitat '{habitat}'. "
            f"Choose from: {', '.join(sorted(VALID_HABITATS))}"
        )

    return None if habitat == "none" else habitat


# --------------------------
# Database lookup
# --------------------------

def habitat_database_names(
    habitat: str,
    config: dict,
) -> List[str]:
    """Return configured Kraken2 database names for a habitat."""
    if habitat not in VALID_HABITATS - {"none"}:
        raise ValueError(f"Unsupported habitat '{habitat}'")

    settings = config.get("Habitat_refinement", {}).get(habitat)

    if settings is None:
        raise ValueError(
            f"Habitat_refinement.{habitat} is not configured"
        )

    databases = settings.get("kraken2", {}).get("databases", {})

    if not databases:
        raise ValueError(
            "No Kraken2 databases configured under "
            f"Habitat_refinement.{habitat}.kraken2.databases"
        )

    return list(databases.keys())


def get_habitat_kraken_database(
    habitat: str,
    database: str,
    config: dict,
) -> str:
    """Return the path of a habitat Kraken2 database."""
    databases = config["Habitat_refinement"][habitat][
        "kraken2"
    ]["databases"]

    if database not in databases:
        raise ValueError(
            f"Unknown Kraken2 database '{database}' "
            f"for habitat '{habitat}'. "
            f"Available databases: {', '.join(databases)}"
        )

    database_path = databases[database].get("path")

    if not database_path:
        raise ValueError(
            f"No database path configured for {habitat}/{database}"
        )

    return database_path


# --------------------------
# Configuration validation
# --------------------------

def validate_habitat_config(config: dict) -> None:
    """Validate the selected habitat and database-specific settings."""
    habitat = selected_habitat(config)

    if habitat is None:
        return

    databases = habitat_database_names(habitat, config)
    settings = config["Habitat_refinement"][habitat]

    for database in databases:
        get_habitat_kraken_database(
            habitat=habitat,
            database=database,
            config=config,
        )

        bracken_settings = settings.get(
            "bracken", {}
        ).get(database)

        if not bracken_settings:
            raise ValueError(
                "No Bracken configuration found under "
                f"Habitat_refinement.{habitat}.bracken.{database}"
            )

        if not bracken_settings.get("level"):
            raise ValueError(
                f"No Bracken level configured for {habitat}/{database}"
            )

        biom_settings = settings.get(
            "kraken_biom", {}
        ).get(database)

        if not biom_settings:
            raise ValueError(
                "No kraken-biom configuration found under "
                f"Habitat_refinement.{habitat}.kraken_biom.{database}"
            )

        if not biom_settings.get("level"):
            raise ValueError(
                "No kraken-biom level configured for "
                f"{habitat}/{database}"
            )


# --------------------------
# Kraken2 settings
# --------------------------

def habitat_kraken_db(
    wildcards: Any,
    *,
    config: dict,
) -> str:
    """Return the database path for a habitat job."""
    return get_habitat_kraken_database(
        habitat=wildcards.habitat,
        database=wildcards.database,
        config=config,
    )


def habitat_kraken_options(
    wildcards: Any,
    *,
    config: dict,
) -> str:
    """Return habitat-specific Kraken2 options."""
    return config["Habitat_refinement"][wildcards.habitat][
        "kraken2"
    ].get("options", "")


def habitat_kraken_threads(
    wildcards: Any,
    *,
    config: dict,
) -> int:
    """Return the requested Kraken2 thread count."""
    return config["Habitat_refinement"][wildcards.habitat][
        "kraken2"
    ].get("threads", 2)


# --------------------------
# Bracken settings
# --------------------------

def habitat_bracken_options(
    wildcards: Any,
    *,
    config: dict,
) -> str:
    """Return habitat-specific Bracken options."""
    return config["Habitat_refinement"][wildcards.habitat][
        "bracken"
    ].get("options", "")


def habitat_bracken_level(
    wildcards: Any,
    *,
    config: dict,
) -> str:
    """Return the Bracken rank option for a habitat database."""
    try:
        return config["Habitat_refinement"][wildcards.habitat][
            "bracken"
        ][wildcards.database]["level"]

    except KeyError as error:
        raise ValueError(
            "No Bracken level configured for "
            f"{wildcards.habitat}/{wildcards.database}"
        ) from error


# --------------------------
# Kraken-BIOM settings
# --------------------------

def habitat_biom_options(
    wildcards: Any,
    *,
    config: dict,
) -> str:
    """Return habitat-specific Kraken-BIOM options."""
    return config["Habitat_refinement"][wildcards.habitat][
        "kraken_biom"
    ].get("options", "")


def habitat_biom_level(
    wildcards: Any,
    *,
    config: dict,
) -> str:
    """Return the Kraken-BIOM rank option for a habitat database."""
    try:
        return config["Habitat_refinement"][wildcards.habitat][
            "kraken_biom"
        ][wildcards.database]["level"]

    except KeyError as error:
        raise ValueError(
            "No kraken-biom level configured for "
            f"{wildcards.habitat}/{wildcards.database}"
        ) from error


# --------------------------
# Bracken report paths
# --------------------------

def habitat_taxonomy_report(
    results_dir: Path,
    sample: str,
    habitat: str,
    database: str,
) -> str:
    """Return one sample's Bracken Kraken-style report."""
    return str(
        results_dir
        / "Habitat_refinement"
        / sample
        / habitat
        / database
        / "bracken"
        / f"{sample}.{database}.bracken.k2report"
    )


def habitat_taxonomy_reports(
    results_dir: Path,
    samples: List[str],
    habitat: str,
    database: str,
) -> List[str]:
    """Return all sample reports for one habitat/database."""
    return [
        habitat_taxonomy_report(
            results_dir=results_dir,
            sample=sample,
            habitat=habitat,
            database=database,
        )
        for sample in samples
    ]


def habitat_reports_for_database(
    wildcards: Any,
    *,
    results_dir: Path,
    samples: List[str],
) -> List[str]:
    """Return the reports required by a habitat summary job."""
    return habitat_taxonomy_reports(
        results_dir=results_dir,
        samples=samples,
        habitat=wildcards.habitat,
        database=wildcards.database,
    )


# --------------------------
# Input reads
# --------------------------

def habitat_read(
    wildcards: Any,
    mate: int,
    *,
    results_dir: Path,
) -> str:
    """Select Soil non-rRNA reads or Glacier glacial-order SSU reads."""
    sample = wildcards.sample

    if wildcards.habitat == "Soil":
        return str(
            results_dir
            / "RNA_Classification"
            / "classified"
            / sample
            / "non_rRNA"
            / f"{sample}_non_rRNA_{mate}.fastq.gz"
        )

    if wildcards.habitat == "Glacier":
        return str(
            results_dir
            / "Habitat_refinement"
            / sample
            / "Glacier"
            / "grouped_reads"
            / f"{sample}.glacial_orders_R{mate}.fastq.gz"
        )

    raise ValueError(
        "No read input mapping defined for habitat "
        f"'{wildcards.habitat}'"
    )


def habitat_read_1(
    wildcards: Any,
    *,
    results_dir: Path,
) -> str:
    """Return the habitat R1 input."""
    return habitat_read(wildcards, 1, results_dir=results_dir)


def habitat_read_2(
    wildcards: Any,
    *,
    results_dir: Path,
) -> str:
    """Return the habitat R2 input."""
    return habitat_read(wildcards, 2, results_dir=results_dir)