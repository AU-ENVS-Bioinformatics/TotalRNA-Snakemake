from pathlib import Path
from typing import List

from snakemake.io import expand


def kraken_database_names(config: dict) -> List[str]:
    databases = config["taxonomy"]["SSU"]["kraken2"].get("databases", {})
    if not databases:
        raise ValueError(
            "No SSU Kraken2 databases configured under "
            "taxonomy -> SSU -> kraken2 -> databases"
        )
    return list(databases.keys())


def get_kraken_database(database: str, config: dict) -> str:
    databases = config["taxonomy"]["SSU"]["kraken2"]["databases"]
    if database not in databases:
        raise ValueError(
            f"Unknown SSU Kraken2 database '{database}'. Available databases: "
            + ", ".join(databases.keys())
        )
    database_path = databases[database].get("path")
    if not database_path:
        raise ValueError(f"No path configured for SSU Kraken2 database '{database}'")
    return database_path


def taxonomy_report(results_dir: Path, sample: str, database: str) -> str:
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
    return [taxonomy_report(results_dir, sample, database) for sample in samples]


def kraken_outputs(
    results_dir: Path,
    samples: List[str],
    config: dict,
) -> List[str]:
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


def bracken_outputs(
    results_dir: Path,
    samples: List[str],
    config: dict,
) -> List[str]:
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


def downstream_outputs(results_dir: Path, config: dict) -> List[str]:
    databases = kraken_database_names(config)
    outputs: List[str] = []
    outputs += expand(
        f"{results_dir}/Taxonomy_Profiling/summary/SSU/{{database}}/taxonomy/{{filename}}",
        database=databases,
        filename=["taxonomy.biom", "taxonomy.tsv"],
    )
    outputs += expand(
        f"{results_dir}/Taxonomy_Profiling/summary/SSU/{{database}}/statistics/{{filename}}",
        database=databases,
        filename=["taxonomy_qualitative.txt", "taxonomy_observations.txt"],
    )
    outputs += expand(
        f"{results_dir}/Taxonomy_Profiling/summary/SSU/{{database}}/phyloseq/{{filename}}",
        database=databases,
        filename=["phyloseq_raw.rds", "phyloseq_filtered.rds"],
    )
    return outputs


def its_taxonomy_enabled(config: dict) -> bool:
    return (
        config["RNA"].get("secondary_method", "none") == "ITS"
        and "ITS" in config.get("taxonomy", {})
    )

def its_blast_regions(config: dict) -> List[str]:
    regions = config["taxonomy"]["ITS"]["blastn"].get(
        "regions", ["full", "ITS1", "5_8S", "ITS2"]
    )
    valid_regions = {"full", "ITS1", "5_8S", "ITS2"}
    invalid_regions = set(regions) - valid_regions
    if invalid_regions:
        raise ValueError(
            "Unsupported ITS BLAST region(s): " + ", ".join(sorted(invalid_regions))
        )
    if len(regions) != len(set(regions)):
        raise ValueError("Duplicate regions configured under taxonomy -> ITS -> blastn -> regions")
    return regions


def its_taxonomy_outputs(
    results_dir: Path,
    samples: List[str],
    config: dict,
) -> List[str]:
    if not its_taxonomy_enabled(config):
        return []
    regions = its_blast_regions(config)
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
        (
            f"{results_dir}/Taxonomy_Profiling/summary/ITS/matrices/"
            f"ITS.assignments.{summary_rank}.counts.tsv"
        ),
        (
            f"{results_dir}/Taxonomy_Profiling/summary/ITS/matrices/"
            "ITS.unassignments.reasons.counts.tsv"
        ),
        (
            f"{results_dir}/Taxonomy_Profiling/summary/ITS/matrices/"
            "ITS.allassignments.status.counts.tsv"
        ),
    ]

    return outputs


def taxonomy_outputs(
    results_dir: Path,
    samples: List[str],
    config: dict,
) -> List[str]:

    # the priority is full > ITS2 > ITS1 > 5_8S
    
    outputs: List[str] = []
    outputs += kraken_outputs(results_dir, samples, config)
    outputs += bracken_outputs(results_dir, samples, config)
    outputs += downstream_outputs(results_dir, config)
    outputs += its_taxonomy_outputs(results_dir, samples, config)
    return outputs
