from typing import List
from pathlib import Path
from snakemake.io import expand


def qc_outputs(qc_dir: Path, samples: List[str]) -> List[str]:
    """
    Return all expected QC output files.

    Parameters
    ----------
    results_dir : Path
        Base results directory.
    samples : List[str]
        List of sample IDs.

    Returns
    -------
    List[str]
        List of expected QC output file paths.
    """

    outputs: List[str] = []

    # FastQC outputs (R1 / R2)
    outputs += expand(
        f"{qc_dir}/{{sample}}/fastqc/{{read}}",
        zip,
        sample=[s for s in samples for _ in ("R1", "R2")],
        read=["R1", "R2"] * len(samples),
    )

    # Fastp reports
    outputs += expand(
        f"{qc_dir}/{{sample}}/trimmed/{{sample}}_fastp.html",
        sample=samples,
    )

    # BWA alignment to contaminants
    outputs += expand(
        f"{qc_dir}/{{sample}}/decontamination/{{sample}}_aligned_contaminants.sam",
        sample=samples,
    )

    # Contaminant read IDs
    outputs += expand(
        f"{qc_dir}/{{sample}}/decontamination/{{sample}}_contaminants_id.txt",
        sample=samples,
    )

    # Cleaned reads (only specifying R1 here, but R2 will be generated in the same step)
    outputs += expand(
        f"{qc_dir}/{{sample}}/decontamination/{{sample}}_R1.cleaned.fastq.gz",
        sample=samples,
    )

    return outputs