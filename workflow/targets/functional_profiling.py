"""Expected completion outputs for functional profiling."""

from pathlib import Path
from typing import List

from snakemake.io import expand

from utils.functional_profiling import (
    SALMON_REFERENCES,
    SALMON_METRICS,
    read_based_enabled,
    assembly_based_enabled,
    validate_functional_config,
)


# --------------------------
# Read-based outputs
# --------------------------

def read_based_functional_profiling_outputs(
    results_dir: Path,
    samples: List[str],
    config: dict,
) -> List[str]:
    """Return MetaPhlAn and optional HUMAnN outputs."""
    if not read_based_enabled(config):
        return []

    functional_dir = results_dir / "Functional_Profiling/read_based"
    merged_dir = functional_dir / "humann_merged"

    outputs: List[str] = []

    outputs += expand(
        f"{functional_dir}/{{sample}}/metaphlan/"
        "{{sample}}_metaphlan_profile.tsv",
        sample=samples,
    )

    if config["functional_profiling"]["read_method"] == "humann":
        outputs += expand(
            f"{functional_dir}/{{sample}}/humann/"
            "{{sample}}_{{metric}}.tsv",
            sample=samples,
            metric=[
                "genefamilies",
                "pathabundance",
                "pathcoverage",
            ],
        )

        outputs += expand(
            f"{merged_dir}/merged_{{metric}}{{suffix}}.tsv",
            metric=[
                "genefamilies",
                "pathabundance",
                "pathcoverage",
            ],
            suffix=[
                "",
                "_stratified",
                "_unstratified",
            ],
        )

        outputs += expand(
            f"{merged_dir}/merged_{{metric}}_relab.tsv",
            metric=["genefamilies", "pathabundance"],
        )

        outputs.append(
            f"{merged_dir}/merged_genefamilies_rxn.tsv"
        )

    return outputs


# --------------------------
# Assembly-based outputs
# --------------------------

def assembly_based_functional_profiling_outputs(
    results_dir: Path,
    samples: List[str],
    config: dict,
) -> List[str]:
    """Return ORF, annotation, and Salmon matrix outputs."""
    if not assembly_based_enabled(config):
        return []

    functional_dir = results_dir / "Functional_Profiling"

    outputs: List[str] = [
        f"{functional_dir}/ORF_prediction/final_orf/ORF_proteins.faa",
        f"{functional_dir}/ORF_prediction/final_orf/ORF_genes.fasta",
        f"{functional_dir}/ORF_prediction/final_orf/ORF_genes.gff3",
    ]

    if config["functional_profiling"]["ORF_predictor"] == "transdecoder":
        outputs += [
            (
                f"{functional_dir}/diamond/transdecoder/"
                "transdecoder.blastp.outfmt6"
            ),
            f"{functional_dir}/hmmscan/merged/pfam.domtblout",
        ]

    outputs.append(
        f"{functional_dir}/eggnog/eggnog.emapper.annotations"
    )

    outputs += expand(
        f"{functional_dir}/salmon/{{reference}}/{{metric}}.tsv",
        reference=SALMON_REFERENCES,
        metric=SALMON_METRICS,
    )

    return outputs


# --------------------------
# Module completion outputs
# --------------------------

def functional_profiling_outputs(
    results_dir: Path,
    samples: List[str],
    config: dict,
) -> List[str]:
    """Return outputs for the selected functional branches."""
    validate_functional_config(config)

    outputs: List[str] = []

    outputs += read_based_functional_profiling_outputs(
        results_dir,
        samples,
        config,
    )

    outputs += assembly_based_functional_profiling_outputs(
        results_dir,
        samples,
        config,
    )

    return outputs