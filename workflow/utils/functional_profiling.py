"""
Helpers for functional profiling.

Provides configuration validation, branch selection, final ORF
source paths, Salmon references, and HMMscan chunk discovery.

Expected completion outputs remain in targets/functional_profiling.py.
"""

from pathlib import Path
from typing import Any, Dict, List


VALID_FUNCTIONAL_MODULES = {
    "read",
    "assembly",
    "both",
    "all",
}

SALMON_REFERENCES = ["transcripts", "cds"]
SALMON_METRICS = ["TPM", "NumReads"]


# --------------------------
# Branch selection
# --------------------------

def read_based_enabled(config: dict) -> bool:
    """Return whether read-based profiling is selected."""
    return config["functional_profiling"]["module"] in {
        "read", "both", "all"
    }


def assembly_based_enabled(config: dict) -> bool:
    """Return whether assembly-based profiling is selected."""
    return config["functional_profiling"]["module"] in {
        "assembly", "both", "all"
    }


# --------------------------
# Configuration validation
# --------------------------

def validate_functional_config(config: dict) -> None:
    """Validate the selected functional profiling branches."""
    settings = config["functional_profiling"]
    module = settings["module"]

    if module not in VALID_FUNCTIONAL_MODULES:
        raise ValueError(
            f"Unsupported functional_profiling.module '{module}'. "
            f"Choose from: {', '.join(sorted(VALID_FUNCTIONAL_MODULES))}"
        )

    if read_based_enabled(config):
        read_method = settings["read_method"]

        if read_method not in {"metaphlan", "humann"}:
            raise ValueError(
                "Unsupported functional_profiling.read_method: "
                f"{read_method}"
            )

    if assembly_based_enabled(config):
        predictor = settings["ORF_predictor"]

        if predictor not in {"prodigal", "transdecoder"}:
            raise ValueError(
                "Unsupported functional_profiling.ORF_predictor: "
                f"{predictor}"
            )


# --------------------------
# Final ORF sources
# --------------------------

def final_orf_sources(
    results_dir: Path,
    config: dict,
) -> Dict[str, str]:
    """Return predictor-specific sources for standardized ORF outputs."""
    predictor = config["functional_profiling"]["ORF_predictor"]
    orf_dir = results_dir / "Functional_Profiling" / "ORF_prediction"

    if predictor == "prodigal":
        return {
            "proteins": str(orf_dir / "prodigal/ORF_proteins.faa"),
            "genes": str(orf_dir / "prodigal/ORF_genes.fasta"),
            "gff": str(orf_dir / "prodigal/ORF_genes.gff3"),
        }

    if predictor == "transdecoder":
        predict_dir = orf_dir / "transdecoder/predict"

        return {
            "proteins": str(
                predict_dir / "transcripts.fasta.transdecoder.pep"
            ),
            "genes": str(
                predict_dir / "transcripts.fasta.transdecoder.cds"
            ),
            "gff": str(
                predict_dir / "transcripts.fasta.transdecoder.gff3"
            ),
        }

    raise ValueError(f"Unsupported ORF predictor: {predictor}")


# --------------------------
# Salmon references
# --------------------------

def salmon_reference(wildcards: Any, *, results_dir: Path) -> str:
    """Return the filtered transcript or standardized CDS reference."""
    if wildcards.reference == "transcripts":
        return str(results_dir / "Assembly/coassembly/non_rRNA/final/transcript.fasta")
    
    if wildcards.reference == "cds":
        return str(results_dir / "Functional_Profiling/ORF_prediction/final_orf/ORF_genes.fasta")

    raise ValueError(f"Unsupported Salmon reference: {wildcards.reference}")

# --------------------------
# HMMscan checkpoint discovery
# --------------------------

def get_hmmscan_outputs(
    *,
    checkpoints: Any,
    results_dir: Path,
) -> List[str]:
    """Return expected domain results for discovered peptide chunks."""
    checkpoint_result = checkpoints.hmmscan_split.get()
    chunk_dir = Path(checkpoint_result.output.chunks)

    chunk_files = sorted(chunk_dir.glob("*.pep"))

    if not chunk_files:
        raise ValueError(
            f"No peptide chunks found in {chunk_dir}"
        )

    return [
        str(
            results_dir
            / "Functional_Profiling"
            / "hmmscan"
            / "results"
            / f"{chunk.stem}.pfam.domtblout"
        )
        for chunk in chunk_files
    ]