from pathlib import Path
from typing import List

from snakemake.io import expand


VALID_FUNCTIONAL_MODULES = {
    "read",
    "assembly",
    "both",
    "all",
}

SALMON_REFERENCES = [
    "transcripts",
    "cds",
]

SALMON_METRICS = [
    "TPM",
    "NumReads",
]


################################################################################
# READ-BASED OUTPUTS
################################################################################

def read_based_functional_profiling_outputs(
    results_dir: Path,
    samples: List[str],
    config: dict,
) -> List[str]:
    """
    Return MetaPhlAn and optional HUMAnN outputs for all samples.
    """

    functional_module = config["functional_profiling"]["module"]

    if functional_module not in ["read", "both", "all"]:
        return []

    read_method = config["functional_profiling"]["read_method"]

    if read_method not in ["metaphlan", "humann"]:
        raise ValueError(
            f"Unsupported functional_profiling.read_method: {read_method}"
        )

    functional_dir = f"{results_dir}/Functional_Profiling/read_based"
    merged_dir = f"{functional_dir}/humann_merged"

    outputs: List[str] = []

    # MetaPhlAn is required for both read methods.
    outputs += expand(
        f"{functional_dir}/{{sample}}/metaphlan/{{sample}}_metaphlan_profile.tsv",
        sample=samples,
    )

    if read_method == "humann":

        # Sample-specific HUMAnN tables.
        outputs += expand(
            f"{functional_dir}/{{sample}}/humann/{{sample}}_{{metric}}.tsv",
            sample=samples,
            metric=[
                "genefamilies",
                "pathabundance",
                "pathcoverage",
            ],
        )

        # Joined tables and their stratified/unstratified versions.
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

        # Relative abundance applies to abundance tables.
        outputs += expand(
            f"{merged_dir}/merged_{{metric}}_relab.tsv",
            metric=[
                "genefamilies",
                "pathabundance",
            ],
        )

        # UniRef90 gene families regrouped into reactions.
        outputs += [
            f"{merged_dir}/merged_genefamilies_rxn.tsv",
        ]

    return outputs


################################################################################
# ASSEMBLY-BASED OUTPUTS
################################################################################

def assembly_based_functional_profiling_outputs(
    results_dir: Path,
    samples: List[str],
    config: dict,
) -> List[str]:
    """
    Return coassembly annotation outputs and Salmon abundance matrices.
    """

    functional_module = config["functional_profiling"]["module"]

    if functional_module not in ["assembly", "both", "all"]:
        return []

    orf_predictor = config["functional_profiling"]["ORF_predictor"]

    if orf_predictor not in ["prodigal", "transdecoder"]:
        raise ValueError(
            f"Unsupported functional_profiling.ORF_predictor: {orf_predictor}"
        )

    functional_dir = f"{results_dir}/Functional_Profiling"

    outputs: List[str] = []

    # Standardized final ORF outputs.
    outputs += [
        f"{functional_dir}/ORF_prediction/final_orf/ORF_proteins.faa",
        f"{functional_dir}/ORF_prediction/final_orf/ORF_genes.fasta",
        f"{functional_dir}/ORF_prediction/final_orf/ORF_genes.gff3",
    ]

    # TransDecoder evidence outputs.
    if orf_predictor == "transdecoder":
        outputs += [
            f"{functional_dir}/diamond/transdecoder/transdecoder.blastp.outfmt6",
            f"{functional_dir}/hmmscan/merged/pfam.domtblout",
        ]

    # eggNOG functional annotation.
    outputs += [
        f"{functional_dir}/eggnog/eggnog.emapper.annotations",
    ]

    # Cross-sample Salmon abundance matrices.
    outputs += expand(
        f"{functional_dir}/salmon/{{reference}}/{{metric}}.tsv",
        reference=SALMON_REFERENCES,
        metric=SALMON_METRICS,
    )

    return outputs


################################################################################
# COMBINED FUNCTIONAL PROFILING TARGETS
################################################################################

def functional_profiling_outputs(
    results_dir: Path,
    samples: List[str],
    config: dict,
) -> List[str]:
    """
    Return outputs for the selected functional profiling branches.

    read       -> read-based outputs
    assembly   -> assembly-based outputs
    both/all   -> both branches
    """

    functional_module = config["functional_profiling"]["module"]

    if functional_module not in VALID_FUNCTIONAL_MODULES:
        raise ValueError(
            f"Unsupported functional_profiling.module '{functional_module}'. "
            f"Choose from: {', '.join(sorted(VALID_FUNCTIONAL_MODULES))}"
        )

    outputs: List[str] = []

    outputs += read_based_functional_profiling_outputs(
        results_dir=results_dir,
        samples=samples,
        config=config,
    )

    outputs += assembly_based_functional_profiling_outputs(
        results_dir=results_dir,
        samples=samples,
        config=config,
    )

    return outputs