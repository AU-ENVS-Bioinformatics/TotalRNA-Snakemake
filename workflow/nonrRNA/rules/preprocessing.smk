##### CONCATENATES READ PAIRS WITHIN ONE SAMPLE TO GENERATE 1 FASTQ

rule concatenate:
    message:
        "[Concatenate] within sample {wildcards.sample} reads"
    input:
        nonrRNA_R1=f"{RESULTS_DIR}/nonrRNA/{{sample}}/filtered/{{sample}}_nonrRNA_1.fastq.gz",
        nonrRNA_R2=f"{RESULTS_DIR}/nonrRNA/{{sample}}/filtered/{{sample}}_nonrRNA_2.fastq.gz",
    output:
        nonrRNA_concatenate=f"{RESULTS_DIR}/nonrRNA/{{sample}}/filtered/{{sample}}_nonrRNA_concat.fastq.gz",
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.nonrRNA_concatenate})

        cat {input.nonrRNA_R1} {input.nonrRNA_R2} > {output.nonrRNA_concatenate}
        """

# consider adding searches for tRNA or similar if proper tools and databases are located
#### CONCATENATE ACROSS SAMPLES
READS = ["1", "2"]

rule concatenate_across_samples:
    input:
        lambda wc: expand(
            f"{RESULTS_DIR}/nonrRNA/{{sample}}/filtered/{{sample}}_nonrRNA_{wc.read}.fastq.gz",
            sample=SAMPLES,
        )
    output:
        merged=f"{RESULTS_DIR}/nonrRNA/concatenated/nonrRNA_{{read}}.fastq.gz"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.merged})

        cat {input} > {output.merged}
        """

"""
In the future consider making it easier to identify read pairs across the various naming conventions

of course this below doesn't fit with the actual output names, but use this as a structure to make something similar
READ1_PATTERNS = [
    "{sample}_1.fastq.gz",
    "{sample}_1.fq.gz",
    "{sample}_1.fastq",
    "{sample}_1.fq",
    "{sample}_R1.fastq.gz",
    "{sample}_R1.fq.gz",
    "{sample}_read1.fastq.gz",
]

READ2_PATTERNS = [
    "{sample}_2.fastq.gz",
    "{sample}_2.fq.gz",
    "{sample}_2.fastq",
    "{sample}_2.fq",
    "{sample}_R2.fastq.gz",
    "{sample}_R2.fq.gz",
    "{sample}_read2.fastq.gz",
]

def find_read(sample, read_patterns):
    for pattern in read_patterns:
        f = Path(RAW_DIR) / pattern.format(sample=sample)
        if f.exists():
            return str(f)

    raise FileNotFoundError(
        f"No matching FASTQ found for sample {sample}"
    )
"""