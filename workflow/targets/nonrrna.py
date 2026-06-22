from typing import List
from pathlib import Path
from snakemake.io import expand

METRICS = ["genefamilies", "pathabundance", "pathcoverage"]

def non_rrna_outputs(
    results_dir: Path,
    samples: List[str],
    config: dict,
) -> List[str]:

    method: str = config["nonrRNA"]["method"]
    outputs: List[str] = []

    if method in ["sequence reads", "metaphlan", "diamond", "humann"]:

        # --------------------------
        # 1. raw per-sample outputs
        # --------------------------
        outputs += expand(
            f"{results_dir}/nonrRNA/{{sample}}/metaphlan/{{sample}}_metaphlan_profile.tsv",
            sample=samples,
        )

        outputs += expand(
            f"{results_dir}/nonrRNA/{{sample}}/humann/{{sample}}_genefamilies.tsv",
            sample=samples,
        )

        # --------------------------
        # 1. link/join input files (IMPORTANT)
        # --------------------------
        outputs += expand(
            f"{results_dir}/nonrRNA/humann_merged/humann_link_input/{{sample}}_{{metric}}.tsv",
            sample=samples,
            metric=METRICS,
        )

        # --------------------------
        # 2. merged (join_tables)
        # --------------------------
        outputs += [
            f"{results_dir}/nonrRNA/humann_merged/merged_{metric}.tsv"
            for metric in METRICS
        ]

        # --------------------------
        # 3. unstratified (MOST IMPORTANT)
        # --------------------------
        outputs += [
            f"{results_dir}/nonrRNA/humann_merged/merged_{metric}_unstratified.tsv"
            for metric in METRICS
        ]

        # --------------------------
        # 4. renormalized
        # --------------------------
        outputs += [
            f"{results_dir}/nonrRNA/humann_merged/merged_{metric}_relab.tsv"
            for metric in METRICS
        ]

        # --------------------------
        # 5. genefamilies KO (optional but common)
        # --------------------------
        outputs += [
            f"{results_dir}/nonrRNA/humann_merged/merged_genefamilies_xrn.tsv"
        ]

        return outputs

    else:
        raise ValueError(f"Unsupported rRNA method: {method}")