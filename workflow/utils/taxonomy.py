"""
Helpers for the taxonomy module.

Provides configuration validation, marker selection, database lookup,
Kraken2/Bracken settings, report paths, and input-read selection.

Configuration and directory dependencies are passed explicitly.
Expected module completion outputs belong in targets/taxonomy.py.
"""

from pathlib import Path
from typing import Any, List

VALID_TAXONOMY_MODULES = {"read_based", "assembly_based", "both"}

READ_BASED_TAXONOMY_MARKERS = {
    "SSU": ("SSU",),
    "ITS": ("ITS",),
    "all": ("SSU", "ITS"),
    "both": ("SSU", "ITS"),
}

def taxonomy_module(config: dict) -> str:
    """Return the selected taxonomy approach."""
    module = config.get("taxonomy", {}).get("module", "read_based")

    if not isinstance(module, str) or module not in VALID_TAXONOMY_MODULES:
        raise ValueError(
            f"Unsupported taxonomy.module {module!r}. "
            "Choose read_based, assembly_based, or both. "
            "Set SSU, ITS, or all under taxonomy.read_based.refinement."
        )

    return module

def read_based_taxonomy_markers(config: dict) -> List[str]:
    """Return the markers requested for read-based taxonomy."""
    if taxonomy_module(config) == "assembly_based":
        return []

    settings = config.get("taxonomy", {}).get("read_based", {})
    if not isinstance(settings, dict):
        raise ValueError("taxonomy.read_based must be a mapping.")

    refinement = settings.get("refinement", "SSU")
    if (
        not isinstance(refinement, str)
        or refinement not in READ_BASED_TAXONOMY_MARKERS
    ):
        raise ValueError(
            f"Unsupported taxonomy.read_based.refinement {refinement!r}. "
            "Choose SSU, ITS, or all."
        )

    return list(READ_BASED_TAXONOMY_MARKERS[refinement])

def taxonomy_markers(config: dict) -> List[str]:
    """Return read-based markers for existing callers."""
    return read_based_taxonomy_markers(config)

def taxonomy_marker_selected(marker: str, config: dict) -> bool:
    return marker in read_based_taxonomy_markers(config)

def validate_taxonomy_config(config: dict) -> None:
    """
    Validate dependencies between taxonomy and upstream modules.
    """

    module = taxonomy_module(config)

    markers = read_based_taxonomy_markers(config)

    for marker in markers:
        if marker not in config["taxonomy"]:
            raise ValueError(
                f"taxonomy.read_based.refinement requests {marker!r}, but taxonomy.{marker} is missing."
            )

    if "ITS" in markers:
        secondary_method = config.get("RNA_classification", {}).get(
            "secondary_method",
            "none",
        )

        if secondary_method != "ITS":
            raise ValueError(
                "ITS taxonomy was selected, but RNA_classification.secondary_method is not 'ITS'"
            )

        reconstruction_targets = config.get(
            "reconstruction",
            {},
        ).get("targets", [])

        if "ITS_candidates" not in reconstruction_targets:
            raise ValueError(
                "ITS taxonomy includes reconstructed ITS BLAST, but 'ITS_candidates' is missing from reconstruction.targets"
            )


#######

def kraken_database_names(
    marker: str,
    config: dict,
) -> List[str]:

    if marker not in {"SSU", "ITS"}:
        raise ValueError(
            f"Unsupported taxonomy marker '{marker}'"
        )

    databases = (
        config["taxonomy"][marker]["kraken2"].get("databases", {})
    )

    # e.g. databases = config["taxonomy"]["SSU"]["kraken2"].get("databases", {})

    if not databases:
        raise ValueError(
            f"No Kraken2 databases configured under taxonomy -> {marker} -> kraken2 -> databases"
        )

    return list(databases.keys())

def get_kraken_database(
    marker: str,
    database: str,
    config: dict,
) -> str:

    databases = (
        config["taxonomy"][marker]["kraken2"]["databases"]
    )

    if database not in databases:
        raise ValueError(
            f"Unknown {marker} Kraken2 database '{database}'. Available databases: {', '.join(databases)}"
        )

    database_path = databases[database].get("path")

    if not database_path:
        raise ValueError(
            f"No path configured for {marker} Kraken2 database '{database}'"
        )

    return database_path


# =========================================================
# Kraken2 wrappers
# =========================================================

def kraken_db(wildcards: Any, *, config: dict) -> str:
    """Return the Kraken2 database path for a job."""
    return get_kraken_database(
        marker=wildcards.marker,
        database=wildcards.database,
        config=config,
    )


def kraken_options(wildcards: Any, *, config: dict) -> str:
    """Return Kraken2 options for the job's marker."""
    return config["taxonomy"][wildcards.marker]["kraken2"].get(
        "options", ""
    )


def kraken_threads(wildcards: Any, *, config: dict) -> int:
    """Return the requested Kraken2 thread count."""
    return config["taxonomy"][wildcards.marker]["kraken2"].get(
        "threads", 2
    )

# --------------------------
# Bracken settings
# --------------------------

def bracken_options(wildcards: Any, *, config: dict) -> str:
    """Return Bracken options for the job's marker."""
    return config["taxonomy"][wildcards.marker]["bracken"].get(
        "options", ""
    )


def bracken_level(wildcards: Any, *, config: dict) -> str:
    """Return the Bracken rank option for the job's database."""
    try:
        return config["taxonomy"][wildcards.marker]["bracken"][
            wildcards.database
        ]["level"]

    except KeyError as error:
        raise ValueError(
            "No Bracken level configured for "
            f"{wildcards.marker}/{wildcards.database}"
        ) from error


# --------------------------
# Kraken-BIOM settings
# --------------------------

def biom_options(wildcards: Any, *, config: dict) -> str:
    """Return Kraken-BIOM options for the job's marker."""
    return config["taxonomy"][wildcards.marker]["kraken_biom"].get(
        "options", ""
    )


def biom_level(wildcards: Any, *, config: dict) -> str:
    """Return the Kraken-BIOM rank option for the database."""
    try:
        return config["taxonomy"][wildcards.marker]["kraken_biom"][
            wildcards.database
        ]["level"]

    except KeyError as error:
        raise ValueError(
            "No kraken-biom level configured for "
            f"{wildcards.marker}/{wildcards.database}"
        ) from error


# --------------------------
# Bracken report paths
# --------------------------

def taxonomy_report(
    results_dir: Path,
    sample: str,
    marker: str,
    database: str,
) -> str:
    """Return a sample's Bracken Kraken-style report path."""
    return str(
        results_dir
        / "Taxonomy_Profiling"
        / "read_based"
        / sample
        / marker
        / database
        / "bracken"
        / f"{sample}.{database}.bracken.k2report"
    )


def taxonomy_reports(
    results_dir: Path,
    samples: List[str],
    marker: str,
    database: str,
) -> List[str]:
    """Return Bracken report paths across samples."""
    return [
        taxonomy_report(
            results_dir=results_dir,
            sample=sample,
            marker=marker,
            database=database,
        )
        for sample in samples
    ]


def taxonomy_reports_for_database(
    wildcards: Any,
    *,
    results_dir: Path,
    samples: List[str],
) -> List[str]:
    """Return reports for the job's marker and database."""
    return taxonomy_reports(
        results_dir=results_dir,
        samples=samples,
        marker=wildcards.marker,
        database=wildcards.database,
    )


# --------------------------
# Taxonomy input reads
# --------------------------

def original_ssu_read(
    wildcards: Any,
    mate: int,
    *,
    results_dir: Path,
) -> str:
    """Return an original SSU read file from RNA classification."""
    return str(
        results_dir
        / "RNA_Classification"
        / "classified"
        / wildcards.sample
        / "SSU"
        / f"{wildcards.sample}_SSU_{mate}.fastq.gz"
    )


def grouped_ssu_read(
    sample: str,
    group: str,
    mate: int,
    *,
    results_dir: Path,
) -> str:
    """Return a read file extracted from a SILVA taxonomic group."""
    return str(
        results_dir
        / "Taxonomy_Profiling"
        / "read_based"
        / sample
        / "SSU"
        / "silva"
        / "grouped_reads"
        / f"{sample}.{group}_R{mate}.fastq.gz"
    )

def assembly_taxonomy_read(
    wildcards: Any,
    mate: int,
    *,
    results_dir: Path,
) -> str:
    """Return the per-sample reads for alignment to the marker co-assembly."""
    if wildcards.marker == "SSU":
        return original_ssu_read(wildcards, mate, results_dir=results_dir)

    raise ValueError(
        f"No co-assembly read input configured for {wildcards.marker!r}."
    )

def taxonomy_read(
    wildcards: Any,
    mate: int,
    *,
    results_dir: Path,
) -> str:
    """Select original SSU, grouped eukaryotic SSU, or ITS reads."""
    if wildcards.marker == "SSU":
        if wildcards.database == "pr2":
            return grouped_ssu_read(
                sample=wildcards.sample,
                group="eukaryotes",
                mate=mate,
                results_dir=results_dir,
            )

        return original_ssu_read(
            wildcards,
            mate,
            results_dir=results_dir,
        )

    if wildcards.marker == "ITS":
        return str(
            results_dir
            / "RNA_Classification"
            / "intermediate"
            / wildcards.sample
            / "ITS_candidates"
            / f"{wildcards.sample}_ITS_candidates_{mate}.fastq.gz"
        )

    raise ValueError(
        f"Unsupported taxonomy marker '{wildcards.marker}'"
    )


def taxonomy_read_1(
    wildcards: Any,
    *,
    results_dir: Path,
) -> str:
    """Return the taxonomy R1 input."""
    return taxonomy_read(wildcards, 1, results_dir=results_dir)


def taxonomy_read_2(
    wildcards: Any,
    *,
    results_dir: Path,
) -> str:
    """Return the taxonomy R2 input."""
    return taxonomy_read(wildcards, 2, results_dir=results_dir)


# --------------------------
# ITS settings
# --------------------------

def its_blast_regions(config: dict) -> List[str]:
    """Return validated ITS BLAST regions."""
    regions = config["taxonomy"]["ITS"]["blastn"].get(
        "regions",
        ["full", "ITS1", "5_8S", "ITS2"],
    )

    valid_regions = {"full", "ITS1", "5_8S", "ITS2"}
    invalid_regions = set(regions) - valid_regions

    if invalid_regions:
        raise ValueError(
            "Unsupported ITS BLAST region(s): "
            + ", ".join(sorted(invalid_regions))
        )

    if len(regions) != len(set(regions)):
        raise ValueError(
            "Duplicate regions configured under "
            "taxonomy -> ITS -> blastn -> regions"
        )

    return regions


def its_summary_rank(config: dict) -> str:
    """Return the validated ITS summary rank."""
    
    rank = config["taxonomy"]["ITS"].get(
        "summary", {}
    ).get("taxonomic_rank", "genus")

    valid_ranks = {
        "kingdom",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species",
    }

    if rank not in valid_ranks:
        raise ValueError(
            f"Unsupported ITS summary rank: {rank}"
        )

    return rank

# =========================================================
# ASSEMBLY BASED TAXONOMY
# =========================================================
def assembly_based_taxonomy_enabled(config: dict) -> bool:
    """Return whether assembly-based taxonomy is requested."""
    return taxonomy_module(config) in {"assembly_based", "both"}

def bwamem_threads(wildcards: Any, *, config: dict) -> int:
    """Return the requested BWA-MEM2 thread count."""
    return config["taxonomy"][wildcards.marker]["bwamem2"].get(
        "threads", 2
    )