from typing import List
from pathlib import Path
from snakemake.io import expand


def reconstruction_outputs(
    results_dir: Path,
    samples: List[str],
    config: dict,
) -> List[str]:
    """
    Return all expected reconstruction output files.
    """

    reconstruction_method = config["reconstruction"]["method"]
    secondary_method = config["RNA"].get(
        "secondary_method",
        "none",
    )

    if reconstruction_method != "rnaspades":
        raise ValueError(
            f"Unsupported reconstruction method: "
            f"{reconstruction_method}"
        )

    outputs: List[str] = []

    # Sample-specific SSU reconstruction
    outputs += expand(
        f"{results_dir}/Assembly/{{sample}}/SSU/rnaspades/transcripts.fasta",
        sample=samples,
    )

    # Optional sample-specific ITS-candidate reconstruction
    if secondary_method == "ITS":
        outputs += expand(
            f"{results_dir}/Assembly/{{sample}}/ITS_candidates/rnaspades/transcripts.fasta",
            sample=samples,
        )

        outputs.append(
            f"{results_dir}/Assembly/concatenated/ITS/cross_sample_ITS.fasta"
        )

    # Cross-sample non-rRNA reconstruction
    outputs.append(
        f"{results_dir}/Assembly/coassembly/non_rRNA/rnaspades/transcripts.fasta"
    )

    return outputs