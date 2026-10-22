rule concatenate_samples_coassembly:
    wildcard_constraints:
        target=CROSS_SAMPLE_TARGET_PATTERN,
        read="1|2"
    input:
        lambda wc: expand(
            f"{RNA_CLASSIFIED_DIR}/{{sample}}/{wc.target}/{{sample}}_{wc.target}_{wc.read}.fastq.gz",
            sample=SAMPLES,
        )
    output:
        merged=(
            f"{ASSEMBLY_DIR}/coassembly/{{target}}/concatenated/{{target}}_{{read}}.fastq.gz"
        )
    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {output.merged})
        
        cat {input} > {output.merged}
        """