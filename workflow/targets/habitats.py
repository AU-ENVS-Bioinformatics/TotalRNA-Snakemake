from pathlib import Path
from typing import List

from snakemake.io import expand

def kraken_database_names(config: dict) -> List[str]:
    databases = config["Habitat_refinement"]["Soil"]["kraken2"].get("databases", {})
    if not databases:
        raise ValueError(
            "No SSU Kraken2 databases configured under taxonomy -> SSU -> kraken2 -> databases"
        )
    return list(databases.keys())

def habitat_taxonomy_reports(
    results_dir: Path,
    samples: List[str],
    habitat: str,
    database: str,
) -> List[str]:
    """Return the habitat-specific Bracken Kraken-format reports."""

    return [
        str(
            Path(results_dir)
            / "Habitat_refinement"
            / sample
            / habitat
            / database
            / "bracken"
            / f"{sample}.{database}.bracken.k2report"
        )
        for sample in samples
    ]

def kraken_outputs_SoilMicroDB(
    results_dir: Path,
    samples: List[str],
    config: dict,
) -> List[str]:

    databases = kraken_database_names(config)

    outputs: List[str] = []

    outputs += expand(
            f"{results_dir}/Habitat_refinement/{{sample}}/Soil/{{database}}/kraken/{{sample}}.{{database}}.k2report",
        sample=samples,
        database=databases,
    )

    outputs += expand(
            f"{results_dir}/Habitat_refinement/{{sample}}/Soil/{{database}}/kraken/{{sample}}.{{database}}.kraken",
        sample=samples,
        database=databases,
    )

    return outputs

def bracken_outputs_SoilMicroDB(
    results_dir: Path,
    samples: List[str],
    config: dict,
) -> List[str]:

    databases = kraken_database_names(config)

    outputs: List[str] = []

    outputs += expand(
            f"{results_dir}/Habitat_refinement/{{sample}}/Soil/{{database}}/bracken/{{sample}}.{{database}}.bracken.k2report",
        sample=samples,
        database=databases,
    )

    outputs += expand(
            f"{results_dir}/Habitat_refinement/{{sample}}/Soil/{{database}}/bracken/{{sample}}.{{database}}.bracken.tsv",
        sample=samples,
        database=databases,
    )

    return outputs

def downstream_outputs_SoilMicroDB(
    results_dir: Path,
    config: dict,
) -> List[str]:

    databases = kraken_database_names(config)

    outputs: List[str] = []

    outputs += expand(f"{results_dir}/Habitat_refinement/summary/Soil/{{database}}/taxonomy/{{filename}}",
        database=databases,
        filename=[
            "taxonomy.biom",
            "taxonomy.tsv",
        ],
    )

    outputs += expand(f"{results_dir}/Habitat_refinement/summary/Soil/{{database}}/phyloseq/{{filename}}",
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

    outputs: List[str] = []

    outputs += kraken_outputs_SoilMicroDB(
        results_dir,
        samples,
        config,
    )

    outputs += bracken_outputs_SoilMicroDB(
        results_dir,
        samples,
        config,
    )

    outputs += downstream_outputs_SoilMicroDB(
        results_dir,
        config,
    )

    return outputs
