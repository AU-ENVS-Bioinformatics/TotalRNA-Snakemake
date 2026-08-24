from typing import List
from pathlib import Path
from snakemake.io import expand
import os
from snakemake.io import glob_wildcards

METRICS = ["genefamilies", "pathabundance", "pathcoverage"]

# translate the sequences to be quantified with salmon into a keyword 
# (used to look up in the dictionary in salmon.smk), to do both original assembled contigs and predicted CDS

SALMON_REFERENCES = {
    "transcripts": "transcripts.fasta",
    "cds": "transcripts.fasta.transdecoder.cds",
}

SALMON_REFERENCE_NAMES = list(SALMON_REFERENCES.keys())

def non_rrna_outputs(
    results_dir: Path,
    samples: List[str],
    config: dict,
) -> List[str]:

    outputs: List[str] = []
    nonrrna_module = config["nonrRNA"]["nonrrna_module"]

    # --------------------------
    # READ-BASED
    # --------------------------
    if nonrrna_module in ["read", "both", "all"]:
        if config["nonrRNA"]["read_method"] in ["metaphlan", "humann"]:
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
    
    if nonrrna_module in ["coassembly", "both", "all"]:

        outputs += [f"{results_dir}/nonrRNA/coassembly/transcripts.fasta"]

        if config["nonrRNA"]["assembly_predictor"] == "transdecoder":
            # DIAMOND Blastp 
            outputs += [f"{results_dir}/nonrRNA/diamond/transdecoder.blastp.outfmt6"]

            # Hmmer scan against Pfam 
            # remember for now hmmscan is run on chunks, so internally it creates temporary files which are then merged into the final output file below
            outputs += [f"{results_dir}/nonrRNA/hmmscan/pfam.domtblout"] 

            # Transdecoder Predict
            outputs += [f"{results_dir}/nonrRNA/predicted/transdecoder.done"]

            # Abundance quantification of samples against assembled contigs and predicted CDS
            outputs += expand(
                f"{results_dir}/nonrRNA/salmon/{{reference}}/{{sample}}/quant.sf",
                sample=samples,
                reference=SALMON_REFERENCE_NAMES,
            )

        #annotations
        outputs += [f"{results_dir}/nonrRNA/eggnog/eggnog.emapper.annotations"]
    
    # --------------------------
    # VALIDATION
    # --------------------------
    if nonrrna_module not in ["read", "coassembly", "both", "all"]:
        raise ValueError(f"Unsupported nonrRNA module: {nonrrna_module}")

    return outputs