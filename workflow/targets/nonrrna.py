from typing import List
from pathlib import Path


SALMON_REFERENCES = [
    "transcripts",
    "cds",
]

SALMON_METRICS = [
    "TPM",
    "NumReads",
]


def non_rrna_outputs(
    results_dir: Path,
    samples: List[str],
    config: dict,
) -> List[str]:
    """
    Return expected assembly-based functional profiling outputs.
    """

    functional_module = config["functional_profiling"]["module"]
    orf_predictor = config["functional_profiling"]["ORF_predictor"]

    if functional_module != "assembly":
        raise ValueError(
            "Only functional_profiling.module: assembly is currently implemented"
        )

    if orf_predictor not in ["prodigal", "transdecoder"]:
        raise ValueError(
            f"Unsupported functional_profiling.ORF_predictor: {orf_predictor}"
        )

    functional_dir = f"{results_dir}/Functional_Profiling"

    outputs: List[str] = []

    # Standardized final ORF outputs
    outputs += [
        f"{functional_dir}/ORF_prediction/final_orf/ORF_proteins.faa",
        f"{functional_dir}/ORF_prediction/final_orf/ORF_genes.fasta",
        f"{functional_dir}/ORF_prediction/final_orf/ORF_genes.gff3",
    ]

    # TransDecoder evidence outputs
    if orf_predictor == "transdecoder":
        outputs += [
            f"{functional_dir}/diamond/transdecoder/transdecoder.blastp.outfmt6",
            f"{functional_dir}/hmmscan/merged/pfam.domtblout",
        ]

    # eggNOG functional annotation
    outputs += [
        f"{functional_dir}/eggnog/eggnog.emapper.annotations",
    ]

    # Salmon abundance matrices
    for reference in SALMON_REFERENCES:
        for metric in SALMON_METRICS:
            outputs.append(
                f"{functional_dir}/salmon/{reference}/{metric}.tsv"
            )

    return outputs