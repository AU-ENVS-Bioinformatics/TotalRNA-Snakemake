"""
Helpers for RNA classification.

Provides configuration validation, SortMeRNA stage handling,
staged input paths, and reference database lookup.

Tool commands remain in the module's rule files.
Expected completion outputs remain in targets/rna_seperation.py.
"""

from pathlib import Path
from typing import Any, Dict, List


# --------------------------
# Configuration validation
# --------------------------

def validate_rna_config(config: dict) -> None:
    """Validate the RNA classification method and refinement."""
    rna_config = config["RNA_classification"]

    method = rna_config["method"]
    refinement = rna_config["refinement"]

    if method == "sortmerna":
        choices = ["staged", "combined"]
        if refinement not in choices:
            raise ValueError(
                f"Unsupported SortMeRNA refinement: {refinement}, choices are: {choices}"
            )

        if refinement == "staged":
            sortmerna_stages(config)

            if rna_config["sortmerna_staged"]["pipeline"] != [
                "ssu",
                "lsu",
            ]:
                raise ValueError(
                    "Staged SortMeRNA postprocessing requires "
                    "pipeline: [ssu, lsu]"
                )

    elif method == "ribodetector":
        choices = ["bbduk"]
        if refinement not in choices:
            raise ValueError(
                f"Unsupported RiboDetector refinement: {refinement}, choices are: {choices}"
            )

    elif method == "bbmap":
        print("BBMap method selected, no refinement required.")
        pass

    else:
        raise ValueError(
            f"Unsupported RNA classification method: {method}"
        )


# --------------------------
# SortMeRNA stage logic
# --------------------------

def sortmerna_stages(config: dict) -> List[str]:
    """
    Return cumulative SortMeRNA stage names.

    Example: ["ssu", "lsu"] becomes ["ssu", "ssu_lsu"].
    """
    pipeline = config["RNA_classification"][
        "sortmerna_staged"
    ]["pipeline"]

    if (
        not isinstance(pipeline, list)
        or not pipeline
        or not all(
            isinstance(stage, str) and stage.strip()
            for stage in pipeline
        )
    ):
        raise ValueError(
            "RNA_classification.sortmerna_staged.pipeline "
            "must be a non-empty list of stage names."
        )

    return [
        "_".join(pipeline[:i + 1])
        for i in range(len(pipeline))
    ]


def previous_stage_map(
    stages: List[str],
) -> Dict[str, str]:
    """Map each stage to the stage producing its input reads."""
    if not stages:
        raise ValueError("SortMeRNA stages cannot be empty.")

    return {
        stages[0]: "decontamination",
        **{
            stages[i]: stages[i - 1]
            for i in range(1, len(stages))
        },
    }


# --------------------------
# SortMeRNA input helpers
# --------------------------

def sortmerna_r1(
    wildcards: Any,
    *,
    results_dir: Path,
    previous_stage: Dict[str, str],
) -> str:
    """Return the R1 input for a staged SortMeRNA job."""
    previous = previous_stage[wildcards.stage]

    if previous == "decontamination":
        return str(
            results_dir
            / "Quality_Control"
            / wildcards.sample
            / "decontamination"
            / f"{wildcards.sample}_R1.cleaned.fastq.gz"
        )

    return str(
        results_dir
        / "RNA_Classification"
        / "intermediate"
        / wildcards.sample
        / "sortmerna"
        / previous
        / f"{wildcards.sample}_{previous}.nonaligned_fwd.fq.gz"
    )


def sortmerna_r2(
    wildcards: Any,
    *,
    results_dir: Path,
    previous_stage: Dict[str, str],
) -> str:
    """Return the R2 input for a staged SortMeRNA job."""
    previous = previous_stage[wildcards.stage]

    if previous == "decontamination":
        return str(
            results_dir
            / "Quality_Control"
            / wildcards.sample
            / "decontamination"
            / f"{wildcards.sample}_R2.cleaned.fastq.gz"
        )

    return str(
        results_dir
        / "RNA_Classification"
        / "intermediate"
        / wildcards.sample
        / "sortmerna"
        / previous
        / f"{wildcards.sample}_{previous}.nonaligned_rev.fq.gz"
    )


# --------------------------
# SortMeRNA database helpers
# --------------------------

def sortmerna_db(
    wildcards: Any,
    *,
    config: dict,
) -> str:
    """Return the reference database for the current stage."""
    database_name = wildcards.stage.split("_")[-1]
    database_key = f"sortmeRNA_{database_name}"

    databases = config["RNA_classification"]["databases"]

    if database_key not in databases:
        raise ValueError(
            "Missing RNA_classification.databases entry: "
            f"{database_key}"
        )

    return databases[database_key]


def sortmerna_db_idx(
    wildcards: Any,
    *,
    config: dict,
) -> str:
    """Return the reference index for the current stage."""
    database_name = wildcards.stage.split("_")[-1]
    database_key = f"sortmeRNA_{database_name}_idx"

    databases = config["RNA_classification"]["databases"]

    if database_key not in databases:
        raise ValueError(
            "Missing RNA_classification.databases entry: "
            f"{database_key}"
        )

    return databases[database_key]


# --------------------------
# Directory helper
# --------------------------

def rrna_dir(
    *paths: str,
    results_dir: Path,
) -> str:
    """Construct a path beneath the existing results/rRNA directory."""
    return str(results_dir.joinpath("rRNA", *paths))