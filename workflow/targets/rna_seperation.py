"""Expected completion outputs for RNA classification."""

from pathlib import Path
from typing import List

from snakemake.io import expand

from utils.rna_classification import validate_rna_config


# --------------------------
# Shared classified outputs
# --------------------------

def filtered_outputs(
    results_dir: Path,
    samples: List[str],
) -> List[str]:
    """Return standardized SSU and non-rRNA paired FASTQ paths."""
    outputs: List[str] = []

    for read in [1, 2]:
        outputs += expand(
            f"{results_dir}/classified/{{sample}}/SSU/{{sample}}_SSU_{read}.fastq.gz",
            sample=samples,
        )

        outputs += expand(
            f"{results_dir}/classified/{{sample}}/non_rRNA/{{sample}}_non_rRNA_{read}.fastq.gz",
            sample=samples,
        )

    return outputs


# --------------------------
# ITS candidate outputs
# --------------------------

def ITS_candidates(
    results_dir: Path,
    samples: List[str],
) -> List[str]:
    """Return paired ITS-candidate FASTQ paths."""
    outputs: List[str] = []

    for read in [1, 2]:
        outputs += expand(
            f"{results_dir}/intermediate/{{sample}}/ITS_candidates/{{sample}}_ITS_candidates_{read}.fastq.gz",
            sample=samples,
        )

    return outputs


# --------------------------
# Module completion outputs
# --------------------------

def rna_outputs(
    results_dir: Path,
    samples: List[str],
    config: dict,
) -> List[str]:
    """
    Return all outputs required to complete RNA classification.

    results_dir is the RNA_Classification directory, not the
    overall results directory.
    """
    validate_rna_config(config)

    rna_config = config["RNA_classification"]

    print(f"[RNA classification] method={rna_config['method']}, refinement={rna_config['refinement']}")

    outputs: List[str] = filtered_outputs(
        results_dir,
        samples,
    )

    if rna_config["secondary_method"] == "ITS":
        outputs += ITS_candidates(
            results_dir,
            samples,
        )

    return outputs