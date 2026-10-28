"""QC input and directory helpers with explicit workflow dependencies."""

from pathlib import Path


def raw_r1(wildcards, *, samplesheet):
    """Return the raw R1 path for a sample."""
    r1 = samplesheet.loc[wildcards.sample, "R1"]
    print(f"[INFO] sample: {wildcards.sample}: read 1 = {r1}")
    return r1


def raw_r2(wildcards, *, samplesheet):
    """Return the raw R2 path for a sample."""
    r2 = samplesheet.loc[wildcards.sample, "R2"]
    print(f"[INFO] sample: {wildcards.sample}: read 2 = {r2}")
    return r2


def sample_qc_dir(wildcards, *, qc_dir):
    """Return the QC output directory for a sample."""
    return Path(qc_dir) / wildcards.sample