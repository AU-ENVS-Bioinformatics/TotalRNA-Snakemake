from __future__ import annotations

from typing import Literal


Decision = Literal["KEEP", "REMOVE"]


def parse_kraken_human_percentage(report_path: str, taxid: str = "9606") -> float:
    """
    Extract percentage of reads assigned to a specific taxonomic ID from a Kraken2 report.

    Parameters
    ----------
    report_path : str
        Path to Kraken2 report file.
    taxid : str
        Taxonomic ID to extract (default: 9606 = Homo sapiens).

    Returns
    -------
    float
        Percentage of reads assigned to the taxon. Returns 0.0 if not found.
    """

    value: float = 0.0

    with open(report_path, "r") as f:
        for line in f:
            cols = line.rstrip().split("\t")

            if len(cols) < 6:
                continue

            if cols[4] == taxid:
                try:
                    value = float(cols[0])
                except ValueError:
                    value = 0.0

    return value


def needs_decontamination(human_pct: float, threshold: float) -> bool:
    """
    Decide whether host decontamination is required.

    Parameters
    ----------
    human_pct : float
        Percentage of reads classified as human.
    threshold : float
        Threshold above which decontamination is triggered.

    Returns
    -------
    bool
        True if decontamination is required.
    """

    return human_pct > threshold


def decide_flag(human_pct: float, threshold: float) -> Decision:
    """
    Convert contamination percentage into a QC decision.

    Parameters
    ----------
    human_pct : float
        Human read percentage.
    threshold : float
        Cutoff threshold.

    Returns
    -------
    Decision
        KEEP or REMOVE.
    """

    return "REMOVE" if needs_decontamination(human_pct, threshold) else "KEEP"


def load_decontam_flag(flag_path: str) -> bool:
    """
    Load Snakemake-generated decontamination flag.

    Parameters
    ----------
    flag_path : str
        Path to flag file.

    Returns
    -------
    bool
        True if REMOVE, False otherwise.
    """

    try:
        with open(flag_path, "r") as f:
            return f.read().strip() == "REMOVE"
    except FileNotFoundError:
        return False