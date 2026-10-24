"""Expected completion outputs for the taxonomy module."""

from pathlib import Path
from typing import List

from snakemake.io import expand

from utils.taxonomy import (
    validate_taxonomy_config,
    taxonomy_markers,
    kraken_database_names,
    its_blast_regions,
    its_summary_rank,
)


# --------------------------
# Kraken2 outputs
# --------------------------

def kraken_outputs(
    results_dir: Path,
    samples: List[str],
    marker: str,
    config: dict,
) -> List[str]:
    """Return Kraken2 reports and read assignments."""
    databases = kraken_database_names(marker, config)
    outputs: List[str] = []

    outputs += expand(
        f"{results_dir}/Taxonomy_Profiling/{{sample}}/{marker}/{{database}}/kraken/{{sample}}.{{database}}.k2report",
        sample=samples,
        database=databases,
    )

    outputs += expand(
        f"{results_dir}/Taxonomy_Profiling/{{sample}}/{marker}/{{database}}/kraken/{{sample}}.{{database}}.kraken",
        sample=samples,
        database=databases,
    )

    return outputs

# --------------------------
# Bracken outputs
# --------------------------

def bracken_outputs(
    results_dir: Path,
    samples: List[str],
    marker: str,
    config: dict,
) -> List[str]:
    """Return Bracken abundance tables and Kraken-style reports."""
    databases = kraken_database_names(marker, config)
    outputs: List[str] = []

    outputs += expand(
        f"{results_dir}/Taxonomy_Profiling/{{sample}}/{marker}/{{database}}/bracken/{{sample}}.{{database}}.bracken.tsv",
        sample=samples,
        database=databases,
    )

    outputs += expand(
        f"{results_dir}/Taxonomy_Profiling/{{sample}}/{marker}/{{database}}/bracken/{{sample}}.{{database}}.bracken.k2report",
        sample=samples,
        database=databases,
    )

    return outputs

# --------------------------
# Combined taxonomy outputs
# --------------------------

def downstream_outputs(
    results_dir: Path,
    marker: str,
    config: dict,
) -> List[str]:
    """Return taxonomy matrices, statistics, and phyloseq objects."""
    databases = kraken_database_names(marker, config)
    outputs: List[str] = []

    outputs += expand(
        f"{results_dir}/Taxonomy_Profiling/summary/{marker}/{{database}}/taxonomy/{{filename}}",
        database=databases,
        filename=[
            "taxonomy.biom",
            "taxonomy.tsv",
        ],
    )

    outputs += expand(
        f"{results_dir}/Taxonomy_Profiling/summary/{marker}/{{database}}/statistics/{{filename}}",
        database=databases,
        filename=[
            "taxonomy_qualitative.txt",
            "taxonomy_observations.txt",
        ],
    )

    outputs += expand(
        f"{results_dir}/Taxonomy_Profiling/summary/{marker}/{{database}}/phyloseq/{{filename}}",
        database=databases,
        filename=["phyloseq.rds"],
    )

    return outputs

# --------------------------
# ITS BLAST outputs
# --------------------------
def its_blast_outputs(
    results_dir: Path,
    samples: List[str],
    config: dict,
) -> List[str]:
    """Return ITS BLAST results, assignments, and summary matrices."""
    regions = its_blast_regions(config)
    summary_rank = its_summary_rank(config)

    outputs: List[str] = []

    outputs += expand(
        f"{results_dir}/Taxonomy_Profiling/{{sample}}/ITS/{{region}}/blastn/{{sample}}.{{region}}.UNITE.blastn.tsv",
        sample=samples,
        region=regions,
    )

    outputs += expand(
        f"{results_dir}/Taxonomy_Profiling/{{sample}}/ITS/{{region}}/taxonomy/{{sample}}.{{region}}.UNITE.assignments.tsv",
        sample=samples,
        region=regions,
    )

    outputs += expand(
        f"{results_dir}/Taxonomy_Profiling/{{sample}}/ITS/combined/{{sample}}.ITS.{{result_type}}.tsv",
        sample=samples,
        result_type=[
            "assignments",
            "unassignments",
            "allasignments",
        ],
    )

    outputs += [
        (
            f"{results_dir}/Taxonomy_Profiling/summary/ITS/matrices/ITS.assignments.{summary_rank}.counts.tsv"
        ),
        (
            f"{results_dir}/Taxonomy_Profiling/summary/ITS/matrices/ITS.unassignments.reasons.counts.tsv"
        ),
        (
            f"{results_dir}/Taxonomy_Profiling/summary/ITS/matrices/ITS.allassignments.status.counts.tsv"
        ),
    ]

    return outputs

# --------------------------
# Module completion outputs
# --------------------------

def taxonomy_outputs(
    results_dir: Path,
    samples: List[str],
    config: dict,
) -> List[str]:
    """
    Return completion outputs for the selected taxonomy markers.

    Read-based taxonomy is requested for every selected marker.
    Assembly-based ITS taxonomy is requested only when
    sample-specific assembly is enabled.
    """
    validate_taxonomy_config(config)

    markers = taxonomy_markers(config)

    assembly_refinement = config["reconstruction"].get(
        "refinement", "coassembly"
    )

    sample_assembly_enabled = assembly_refinement in {
        "samples", "both", "all"
    }

    outputs: List[str] = []

    # Read-based taxonomy: independent of assembly refinement.
    for marker in markers:
        outputs += kraken_outputs(
            results_dir=results_dir,
            samples=samples,
            marker=marker,
            config=config,
        )

        outputs += bracken_outputs(
            results_dir=results_dir,
            samples=samples,
            marker=marker,
            config=config,
        )

        outputs += downstream_outputs(
            results_dir=results_dir,
            marker=marker,
            config=config,
        )

    # Assembly-based ITS taxonomy requires sample assemblies.
    if "ITS" in markers and sample_assembly_enabled:
        outputs += its_blast_outputs(
            results_dir=results_dir,
            samples=samples,
            config=config,
        )

    return outputs