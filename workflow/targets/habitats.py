from pathlib import Path
from typing import List, Optional

from snakemake.io import expand


VALID_HABITATS = {
    "Soil",
    "Glacier",
    "none",
}


def selected_habitat(config: dict) -> Optional[str]:
    """
    Return the habitat selected under Habitat_refinement.habitat.

    Soil    -> run the Soil branch
    Glacier -> run the Glacier branch
    none    -> disable habitat refinement
    """

    habitat = (
        config.get("Habitat_refinement", {}).get("habitat", "none")
    )

    if habitat not in VALID_HABITATS:
        raise ValueError(
            f"Unsupported Habitat_refinement.habitat '{habitat}'. Choose from: {', '.join(sorted(VALID_HABITATS))}"
        )

    if habitat == "none":
        return None

    return habitat


def habitat_database_names(
    habitat: str,
    config: dict,
) -> List[str]:
    """
    Return all Kraken2 database names configured for a habitat.
    """

    if habitat not in VALID_HABITATS - {"none"}:
        raise ValueError(
            f"Unsupported habitat '{habitat}'"
        )

    habitat_settings = config.get("Habitat_refinement", {}).get(habitat)

    if habitat_settings is None:
        raise ValueError(
            f"Habitat_refinement.habitat requests '{habitat}', but Habitat_refinement.{habitat} is not configured"
        )

    databases = habitat_settings.get("kraken2", {}).get("databases", {})

    if not databases:
        raise ValueError(
            f"No Kraken2 databases configured under Habitat_refinement -> {habitat} -> kraken2 -> databases"
        )

    return list(databases.keys())

def get_habitat_kraken_database(
    habitat: str,
    database: str,
    config: dict,
) -> str:
    """
    Return the path belonging to one habitat Kraken2 database.
    """

    databases = config["Habitat_refinement"][habitat]["kraken2"]["databases"]

    if database not in databases:
        raise ValueError(
            f"Unknown Kraken2 database '{database}' for habitat '{habitat}'. Available databases: {', '.join(databases)}"
        )

    database_path = databases[database].get("path")

    if not database_path:
        raise ValueError(
            f"No path configured for Kraken2 database '{database}' under habitat '{habitat}'"
        )

    return database_path

def validate_habitat_config(config: dict) -> None:
    """
    Validate the selected habitat and its database-specific
    Kraken2, Bracken and kraken-biom settings.
    """

    habitat = selected_habitat(config)

    if habitat is None:
        return

    databases = habitat_database_names(
        habitat=habitat,
        config=config,
    )

    habitat_settings = (
        config["Habitat_refinement"][habitat]
    )

    for database in databases:
        get_habitat_kraken_database(
            habitat=habitat,
            database=database,
            config=config,
        )

        bracken_settings = habitat_settings.get("bracken", {}).get(database)

        if bracken_settings is None:
            raise ValueError(
                f"No Bracken configuration found under Habitat_refinement.{habitat}.bracken.{database}"
            )

        if not bracken_settings.get("level"):
            raise ValueError(
                f"No Bracken level configured for {habitat}/{database}"
            )

        biom_settings = habitat_settings.get("kraken_biom", {}).get(database)

        if biom_settings is None:
            raise ValueError(
                f"No kraken-biom configuration found under Habitat_refinement.{habitat}.kraken_biom.{database}"
            )

        if not biom_settings.get("level"):
            raise ValueError(
                f"No kraken-biom level configured for {habitat}/{database}"
            )


def habitat_taxonomy_report(
    results_dir: Path,
    sample: str,
    habitat: str,
    database: str,
) -> str:
    """
    Return one sample-specific Bracken kreport path.
    """

    return str(
        Path(results_dir)
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
    """
    Return all Bracken kreports for one habitat/database.
    """

    return [
        habitat_taxonomy_report(
            results_dir=results_dir,
            sample=sample,
            habitat=habitat,
            database=database,
        )
        for sample in samples
    ]


def habitat_kraken_outputs(
    results_dir: Path,
    samples: List[str],
    habitat: str,
    config: dict,
) -> List[str]:

    databases = habitat_database_names(
        habitat=habitat,
        config=config,
    )

    outputs: List[str] = []

    outputs += expand(
        f"{results_dir}/Habitat_refinement/"
        f"{{sample}}/{habitat}/{{database}}/kraken/"
        f"{{sample}}.{{database}}.k2report",
        sample=samples,
        database=databases,
    )

    outputs += expand(
        f"{results_dir}/Habitat_refinement/"
        f"{{sample}}/{habitat}/{{database}}/kraken/"
        f"{{sample}}.{{database}}.kraken",
        sample=samples,
        database=databases,
    )

    return outputs


def habitat_bracken_outputs(
    results_dir: Path,
    samples: List[str],
    habitat: str,
    config: dict,
) -> List[str]:

    databases = habitat_database_names(
        habitat=habitat,
        config=config,
    )

    outputs: List[str] = []

    outputs += expand(
        f"{results_dir}/Habitat_refinement/"
        f"{{sample}}/{habitat}/{{database}}/bracken/"
        f"{{sample}}.{{database}}.bracken.tsv",
        sample=samples,
        database=databases,
    )

    outputs += expand(
        f"{results_dir}/Habitat_refinement/"
        f"{{sample}}/{habitat}/{{database}}/bracken/"
        f"{{sample}}.{{database}}.bracken.k2report",
        sample=samples,
        database=databases,
    )

    return outputs


def habitat_downstream_outputs(
    results_dir: Path,
    habitat: str,
    config: dict,
) -> List[str]:

    databases = habitat_database_names(
        habitat=habitat,
        config=config,
    )

    outputs: List[str] = []

    outputs += expand(
        f"{results_dir}/Habitat_refinement/"
        f"summary/{habitat}/{{database}}/taxonomy/"
        f"{{filename}}",
        database=databases,
        filename=[
            "taxonomy.biom",
            "taxonomy.tsv",
        ],
    )

    outputs += expand(
        f"{results_dir}/Habitat_refinement/"
        f"summary/{habitat}/{{database}}/phyloseq/"
        f"{{filename}}",
        database=databases,
        filename=[
            "phyloseq.rds",
        ],
    )

    return outputs


def habitat_refinement_outputs(
    results_dir: Path,
    samples: List[str],
    config: dict,
) -> List[str]:
    """
    Return all requested habitat-refinement outputs.

    If habitat is 'none', return an empty list. The completion
    rule can consequently create its flag without running
    Kraken2, Bracken or downstream processing.
    """

    validate_habitat_config(config)

    habitat = selected_habitat(config)

    if habitat is None:
        return []

    outputs: List[str] = []

    outputs += habitat_kraken_outputs(
        results_dir=results_dir,
        samples=samples,
        habitat=habitat,
        config=config,
    )

    outputs += habitat_bracken_outputs(
        results_dir=results_dir,
        samples=samples,
        habitat=habitat,
        config=config,
    )

    outputs += habitat_downstream_outputs(
        results_dir=results_dir,
        habitat=habitat,
        config=config,
    )

    return outputs