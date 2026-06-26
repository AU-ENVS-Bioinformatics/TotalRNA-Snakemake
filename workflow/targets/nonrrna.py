from typing import List
from pathlib import Path
from snakemake.io import expand

METRICS = ["genefamilies", "pathabundance", "pathcoverage"]

def non_rrna_outputs(
    results_dir: Path,
    samples: List[str],
    config: dict,
) -> List[str]:

    outputs: List[str] = []

    if config["nonrRNA"]["module"] == "read":
        if config["nonrRNA"]["method"] in ["metaphlan", "humann"]:
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
    elif config["nonrRNA"]["module"] == "coassembly":

        outputs += [f"{results_dir}/nonrRNA/coassembly/transcripts.fasta"]

        if config["nonrRNA"]["predictor"] == "transdecoder":
            # DIAMOND Blastp 
            outputs += [f"{results_dir}/nonrRNA/diamond/transdecoder.blastp.outfmt6"]

            # Hmmer scan against Pfam 
            outputs += [f"{results_dir}/nonrRNA/hmmscan/pfam.domtblout"]

            # Transdecoder Predict
            outputs += [f"{results_dir}/nonrRNA/predicted/transdecoder.done"]
        
        #annotations
        outputs += [f"{results_dir}/nonrRNA/eggnog/eggnog.emapper.annotations"]

        return outputs
    else:
        raise ValueError(f"Unsupported rRNA module")
    