rule concatenate:
    message:
        "[metaphlan] species-level microbial profiling for {wildcards.sample}"
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