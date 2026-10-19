"""Expected completion outputs for habitat-specific taxonomy."""

from pathlib import Path
from typing import List

from snakemake.io import expand

from utils.habitats import (
    selected_habitat,
    habitat_database_names,
    validate_habitat_config,
)


def habitat_kraken_outputs(
    results_dir: Path,
    samples: List[str],
    habitat: str,
    config: dict,
) -> List[str]:
    """Return Kraken2 reports and read assignments."""
    databases = habitat_database_names(habitat, config)
    outputs: List[str] = []

    outputs += expand(
        f"{results_dir}/Habitat_refinement/{{sample}}/{habitat}/{{database}}/kraken/{{sample}}.{{database}}.k2report",
        sample=samples,
        database=databases,
    )

    outputs += expand(
        f"{results_dir}/Habitat_refinement/{{sample}}/{habitat}/{{database}}/kraken/{{sample}}.{{database}}.kraken",
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
    """Return Bracken abundance tables and Kraken-style reports."""
    databases = habitat_database_names(habitat, config)
    outputs: List[str] = []

    outputs += expand(
        f"{results_dir}/Habitat_refinement/{{sample}}/{habitat}/{{database}}/bracken/{{sample}}.{{database}}.bracken.tsv",
        sample=samples,
        database=databases,
    )

    outputs += expand(
        f"{results_dir}/Habitat_refinement/{{sample}}/{habitat}/{{database}}/bracken/{{sample}}.{{database}}.bracken.k2report",
        sample=samples,
        database=databases,
    )

    return outputs


def habitat_downstream_outputs(
    results_dir: Path,
    habitat: str,
    config: dict,
) -> List[str]:
    """Return combined taxonomy tables and phyloseq objects."""
    databases = habitat_database_names(habitat, config)
    outputs: List[str] = []

    outputs += expand(
        f"{results_dir}/Habitat_refinement/summary/{habitat}/{{database}}/taxonomy/{{filename}}",
        database=databases,
        filename=["taxonomy.biom", "taxonomy.tsv"],
    )

    outputs += expand(
        f"{results_dir}/Habitat_refinement/summary/{habitat}/{{database}}/phyloseq/phyloseq.rds",
        database=databases,
    )

    return outputs


def habitat_refinement_outputs(
    results_dir: Path,
    samples: List[str],
    config: dict,
) -> List[str]:
    """Return requested outputs, or an empty list when disabled."""
    validate_habitat_config(config)

    habitat = selected_habitat(config)

    if habitat is None:
        return []

    outputs: List[str] = []

    outputs += habitat_kraken_outputs(
        results_dir, samples, habitat, config
    )

    outputs += habitat_bracken_outputs(
        results_dir, samples, habitat, config
    )

    outputs += habitat_downstream_outputs(
        results_dir, habitat, config
    )

    return outputs