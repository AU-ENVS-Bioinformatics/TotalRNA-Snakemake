"""
Helpers for the reconstruction module.

Provides:
    - Reconstruction configuration validation.
    - Sample-specific assembly and concatenation targets.
    - Barrnap kingdom selection.
    - Input routing for sample-specific assemblies.
    - Input routing for cross-sample FASTA concatenation.
    - Parsing of the Barrnap kingdom status report.

Expected completion outputs remain in targets/reconstruction.py.
Rules and checkpoint calls remain in Snakemake files.
"""

import csv
from pathlib import Path
from typing import Any, List


# --------------------------
# Co-assembly routing
# --------------------------

# Add methods here only once their rule files are implemented.
IMPLEMENTED_COASSEMBLY_METHODS = {
    "non_rRNA": {"rnaspades"},
    "SSU": {"rnaspades", "metarib"},
}

COASSEMBLY_METHOD_KEYS = {
    "non_rRNA": "non_rRNA_method",
    "SSU": "rRNA_method",
}


def coassembly_targets(config: dict) -> List[str]:
    """Return the configured cross-sample assembly targets."""
    targets = config["reconstruction"].get(
        "coassembly_targets",
        ["non_rRNA"],
    )

    if (
        not isinstance(targets, list)
        or not targets
        or any(
            not isinstance(target, str)
            or target not in COASSEMBLY_METHOD_KEYS
            for target in targets
        )
    ):
        raise ValueError(
            "reconstruction.coassembly_targets must be a non-empty "
            "list containing non_rRNA and/or SSU."
        )

    if len(targets) != len(set(targets)):
        raise ValueError(
            "reconstruction.coassembly_targets contains duplicates."
        )

    # non-rRNA is mandatory for the requested co-assembly workflow.
    if "non_rRNA" not in targets:
        raise ValueError(
            "reconstruction.coassembly_targets must include non_rRNA."
        )

    return targets


def coassembly_method(config: dict, target: str) -> str:
    """Return the selected assembler for a co-assembly target."""
    if target not in COASSEMBLY_METHOD_KEYS:
        raise ValueError(
            f"Unsupported co-assembly target: {target}"
        )

    key = COASSEMBLY_METHOD_KEYS[target]
    method = config["reconstruction"].get(key)

    supported = IMPLEMENTED_COASSEMBLY_METHODS[target]

    if not isinstance(method, str) or method not in supported:
        raise ValueError(
            f"Unsupported reconstruction.{key}: {method!r}. "
            f"Implemented methods for {target}: "
            f"{', '.join(sorted(supported))}."
        )

    return method


def coassembly_fasta(
    *,
    results_dir: Path,
    config: dict,
    target: str,
) -> str:
    """Return the selected assembler's final co-assembly FASTA."""
    method = coassembly_method(config, target)

    return str(
        results_dir
        / "Assembly"
        / "coassembly"
        / target
        / method
        / "transcripts.fasta"
    )


# --------------------------
# Define postprocessing output to fit both co-assembly and cross-sample SUS
# --------------------------

def assembly_postprocessing_targets(config: dict) -> List[str]:
    """Return targets needed by either assembly workflow."""
    refinement = config["reconstruction"].get("refinement", "coassembly")
    targets = []

    if refinement in {"coassembly", "both", "all"}:
        targets.extend(coassembly_targets(config))

    if refinement in {"samples", "both", "all"} and "SSU" not in targets:
        targets.append("SSU")

    return targets


def assembly_postprocessing_input(wildcards: Any, *, results_dir: Path, config: dict) -> str:
    """Route raw co-assemblies or concatenated sample-specific SSU."""
    if wildcards.assembly == "coassembly":
        if wildcards.target not in coassembly_targets(config):
            raise ValueError(f"Co-assembly target is not enabled: {wildcards.target}")

        return str(results_dir / "Assembly" / "coassembly" / wildcards.target / f"{wildcards.target}_transcripts.fasta")

    if wildcards.assembly == "cross_assembly" and wildcards.target == "SSU":
        return str(results_dir / "Assembly" / "cross_assembly" / "concatenated" / "SSU" / "cross_sample_SSU.fasta")

    raise ValueError(f"Unsupported postprocessing branch: {wildcards.assembly}/{wildcards.target}")

# --------------------------
# Configuration validation
# --------------------------

def validate_reconstruction_config(config: dict) -> None:
    """Validate reconstruction configuration."""
    reconstruction = config["reconstruction"]

    refinement = reconstruction.get("refinement", "coassembly")

    if refinement not in {"coassembly", "samples", "both", "all"}:
        raise ValueError(
            "reconstruction.refinement must be "
            "coassembly, samples, both, or all."
        )

    for target in coassembly_targets(config):
        coassembly_method(config, target)

    secondary_method = config["RNA_classification"].get(
        "secondary_method", "none"
    )

    if secondary_method not in {"none", "ITS"}:
        raise ValueError(
            "Unsupported RNA_classification.secondary_method: "
            f"{secondary_method}. Supported values: none, ITS"
        )


# --------------------------
# Reconstruction targets
# --------------------------

def sample_reconstruction_targets(config: dict) -> List[str]:
    """Return the sample-specific assembly branches."""
    targets: List[str] = ["SSU"]

    if config["RNA_classification"].get(
        "secondary_method", "none"
    ) == "ITS":
        targets.append("ITS_candidates")

    return targets


def concatenation_targets(config: dict) -> List[str]:
    """Return the cross-sample FASTA concatenation branches."""
    targets: List[str] = ["SSU"]

    if config["RNA_classification"].get(
        "secondary_method", "none"
    ) == "ITS":
        targets.append("ITS")

    return targets

# --------------------------
# Sample-specific input reads
# --------------------------

def rnaspades_sample_input(
    wildcards: Any,
    read: int,
    *,
    results_dir: Path,
) -> str:
    """Return the SSU or ITS-candidate FASTQ for an assembly."""
    sample = wildcards.sample
    target = wildcards.target

    if target == "SSU":
        return str(
            results_dir
            / "RNA_Classification"
            / "classified"
            / sample
            / "SSU"
            / f"{sample}_SSU_{read}.fastq.gz"
        )

    if target == "ITS_candidates":
        return str(
            results_dir
            / "RNA_Classification"
            / "intermediate"
            / sample
            / "ITS_candidates"
            / f"{sample}_ITS_candidates_{read}.fastq.gz"
        )

    raise ValueError(
        f"Unsupported sample reconstruction target: {target}"
    )


# --------------------------
# Cross-sample FASTA inputs
# --------------------------

def prefixed_reconstructed_fastas(wildcards: Any, *, results_dir: Path, samples: List[str]) -> List[str]:
    """Return prefixed sample FASTAs for SSU or ITS concatenation."""
    assembly_dir = results_dir / "Assembly" / "cross_assembly"

    if wildcards.target == "SSU":
        return [
            str(assembly_dir / sample / "SSU" / "rnaspades" / "transcripts_prefixed.fasta")
            for sample in samples
        ]

    if wildcards.target == "ITS":
        return [
            str(assembly_dir / sample / "ITS_candidates" / "ITSx" / f"{sample}_ITSx.full_prefixed.fasta")
            for sample in samples
        ]

    raise ValueError(f"Unsupported concatenation target: {wildcards.target}")
