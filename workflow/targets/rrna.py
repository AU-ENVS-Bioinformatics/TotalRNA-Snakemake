from typing import List
from pathlib import Path
from snakemake.io import expand

def rrna_outputs(
    results_dir: Path,
    samples: List[str],
    config: dict,
) -> List[str]:
    """
    Return expected rRNA classification outputs.
    """

    method: str = config["rRNA"]["method"]
    outputs: List[str] = []

    # --------------------------
    # Kraken / Bracken
    # --------------------------
    if method in ["kraken2", "bracken", "taxonomic_composition"]:

        # Per-sample kraken outputs
        outputs += expand(
            f"{results_dir}/rRNA/{{sample}}/classification/{{sample}}.k2report",
            sample=samples,
        )

        # Per-sample bracken outputs
        outputs += expand(
            f"{results_dir}/rRNA/{{sample}}/classification/{{sample}}.bracken.k2report",
            sample=samples,
        )

        # Global aggregated outputs
        outputs += [
            f"{results_dir}/rRNA/taxonomy/bracken_species.biom",
            f"{results_dir}/rRNA/taxonomy/bracken_species_qualitative.txt",
            f"{results_dir}/rRNA/taxonomy/bracken_species_observations.txt",
            f"{results_dir}/rRNA/taxonomy/bracken_species.tsv",
            f"{results_dir}/rRNA/taxonomy/phyloseq.rds",
        ]

        return outputs

    else:
        raise ValueError(f"Unsupported rRNA method: {method}")