from typing import List
from pathlib import Path
from snakemake.io import expand

def non_rrna_outputs(
    results_dir: Path,
    samples: List[str],
    config: dict,
) -> List[str]:
    """
    Return expected non-rRNA functional annotation and profiling
    """

    method: str = config["nonrRNA"]["method"]
    outputs: List[str] = []

    # --------------------------
    # Sequence read based
    # --------------------------
    if method in ["sequence reads", "metaphlan", "diamond", "humann"]:

        # Per-sample metaphlan species-level microbial profiling
        outputs += expand(
            f"{results_dir}/nonrRNA/{{sample}}/metaphlan/{{sample}}_1_metaphlan_profile.tsv",
            sample=samples,
        )

        # Per-sample humann functional profiling
        outputs += expand(
            f"{results_dir}/nonrRNA/{{sample}}/humann/{{sample}}_2_genefamilies.tsv",
            sample=samples,
        )

        return outputs

    else:
        raise ValueError(f"Unsupported rRNA method: {method}")