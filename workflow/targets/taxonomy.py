from pathlib import Path
from typing import List

from snakemake.io import expand

VALID_TAXONOMY_MODULES = {
    "SSU",
    "ITS",
    "both",
    "all",
}

def taxonomy_markers(config: dict) -> List[str]:
    """
    Convert taxonomy.module into the marker branches that should run.

    SSU        -> SSU only
    ITS        -> ITS only
    both/all   -> SSU and ITS
    """

    module = config.get("taxonomy", {}).get("module", "SSU")

    if module not in VALID_TAXONOMY_MODULES:
        raise ValueError(
            f"Unsupported taxonomy.module '{module}'. "
            f"Choose from: {', '.join(sorted(VALID_TAXONOMY_MODULES))}"
        )

    if module == "SSU":
        return ["SSU"]

    if module == "ITS":
        return ["ITS"]

    return ["SSU", "ITS"]


def taxonomy_marker_selected(marker: str, config: dict) -> bool:
    return marker in taxonomy_markers(config)


def validate_taxonomy_config(config: dict) -> None:
    """
    Validate dependencies between taxonomy and upstream modules.
    """

    markers = taxonomy_markers(config)

    for marker in markers:
        if marker not in config["taxonomy"]:
            raise ValueError(
                f"taxonomy.module requests '{marker}', but taxonomy.{marker} is not configured"
            )

    if "ITS" in markers:
        secondary_method = config.get("RNA", {}).get(
            "secondary_method",
            "none",
        )

        if secondary_method != "ITS":
            raise ValueError(
                "ITS taxonomy was selected, but RNA.secondary_method is not 'ITS'"
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

def taxonomy_report(
    results_dir: Path,
    sample: str,
    marker: str,
    database: str,
) -> str:

    return str(
        Path(results_dir)
        / "Taxonomy_Profiling"
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

    return [
        taxonomy_report(
            results_dir=results_dir,
            sample=sample,
            marker=marker,
            database=database,
        )
        for sample in samples
    ]

def kraken_outputs(
    results_dir: Path,
    samples: List[str],
    marker: str,
    config: dict,
) -> List[str]:

    databases = kraken_database_names(marker, config)

    outputs: List[str] = []

    outputs += expand(f"{results_dir}/Taxonomy_Profiling/{{sample}}/{marker}/{{database}}/kraken/{{sample}}.{{database}}.k2report",
        sample=samples,
        database=databases,
    )

    outputs += expand(f"{results_dir}/Taxonomy_Profiling/{{sample}}/{marker}/{{database}}/kraken/{{sample}}.{{database}}.kraken",
        sample=samples,
        database=databases,
    )

    return outputs

def bracken_outputs(
    results_dir: Path,
    samples: List[str],
    marker: str,
    config: dict,
) -> List[str]:

    databases = kraken_database_names(marker, config)

    outputs: List[str] = []

    outputs += expand(f"{results_dir}/Taxonomy_Profiling/{{sample}}/{marker}/{{database}}/bracken/{{sample}}.{{database}}.bracken.tsv",
        sample=samples,
        database=databases,
    )

    outputs += expand(f"{results_dir}/Taxonomy_Profiling/{{sample}}/{marker}/{{database}}/bracken/{{sample}}.{{database}}.bracken.k2report",
        sample=samples,
        database=databases,
    )

    return outputs

def downstream_outputs(
    results_dir: Path,
    marker: str,
    config: dict,
) -> List[str]:

    databases = kraken_database_names(marker, config)

    outputs: List[str] = []

    outputs += expand(f"{results_dir}/Taxonomy_Profiling/summary/{marker}/{{database}}/taxonomy/{{filename}}",
        database=databases,
        filename=[
            "taxonomy.biom",
            "taxonomy.tsv",
        ],
    )

    outputs += expand(f"{results_dir}/Taxonomy_Profiling/summary/{marker}/{{database}}/statistics/{{filename}}",
        database=databases,
        filename=[
            "taxonomy_qualitative.txt",
            "taxonomy_observations.txt",
        ],
    )

    outputs += expand(f"{results_dir}/Taxonomy_Profiling/summary/{marker}/{{database}}/phyloseq/{{filename}}",
        database=databases,
        filename=[
            "phyloseq.rds",
        ],
    )
        #"phyloseq_filtered.rds",

    return outputs

def its_blast_regions(config: dict) -> List[str]:

    regions = (
        config["taxonomy"]["ITS"]["blastn"]
        .get(
            "regions",
            ["full", "ITS1", "5_8S", "ITS2"],
        )
    )

    valid_regions = {
        "full",
        "ITS1",
        "5_8S",
        "ITS2",
    }

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


def its_blast_outputs(
    results_dir: Path,
    samples: List[str],
    config: dict,
) -> List[str]:

    regions = its_blast_regions(config)

    outputs: List[str] = []

    outputs += expand(f"{results_dir}/Taxonomy_Profiling/{{sample}}/ITS/{{region}}/blastn/{{sample}}.{{region}}.UNITE.blastn.tsv",
        sample=samples,
        region=regions,
    )

    outputs += expand(f"{results_dir}/Taxonomy_Profiling/{{sample}}/ITS/{{region}}/taxonomy/{{sample}}.{{region}}.UNITE.assignments.tsv",
        sample=samples,
        region=regions,
    )

    outputs += expand(f"{results_dir}/Taxonomy_Profiling/{{sample}}/ITS/combined/{{sample}}.ITS.{{result_type}}.tsv",
        sample=samples,
        result_type=[
            "assignments",
            "unassignments",
            "allasignments",
        ],
    )

    summary_rank = (
        config["taxonomy"]["ITS"]
        .get("summary", {})
        .get("taxonomic_rank", "genus")
    )

    valid_ranks = {
        "kingdom",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species",
    }

    if summary_rank not in valid_ranks:
        raise ValueError(
            f"Unsupported ITS summary rank: {summary_rank}"
        )

    outputs += [
        f"{results_dir}/Taxonomy_Profiling/summary/ITS/matrices/ITS.assignments.{summary_rank}.counts.tsv",
        f"{results_dir}/Taxonomy_Profiling/summary/ITS/matrices/ITS.unassignments.reasons.counts.tsv",
        f"{results_dir}/Taxonomy_Profiling/summary/ITS/matrices/ITS.allassignments.status.counts.tsv",
    ]

    return outputs


def taxonomy_outputs(
    results_dir: Path,
    samples: List[str],
    config: dict,
) -> List[str]:

    validate_taxonomy_config(config)

    markers = taxonomy_markers(config)
    outputs: List[str] = []

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

    if "ITS" in markers:
        outputs += its_blast_outputs(
            results_dir=results_dir,
            samples=samples,
            config=config,
        )

    return outputs