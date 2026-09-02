READS = ["1", "2"]

rule concatenate_across_samples2:
    input:
        lambda wc: expand(
            f"{RNA_CLASSIFIED_DIR}/{{sample}}/non_rRNA/{{sample}}_non_rRNA_{wc.read}.fastq.gz",
            sample=SAMPLES,
        )
    output:
        merged=f"{ASSEMBLY_DIR}/coassembly/non_rRNA/concatenated/nonrRNA_{{read}}.fastq.gz"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.merged})

        cat {input} > {output.merged}
        """
