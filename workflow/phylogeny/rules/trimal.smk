################################################################################
# 8. Trim alignment with trimAl
################################################################################

rule trimal_ssu:
    conda:
        "../envs/trimal.yaml"
    message:
        "[trimAl] trim SSU alignment for {wildcards.sample}"
    input:
        alignment=f"{PHYLOGENY_DIR}/{{sample}}/Alignment/{{sample}}_SSU_mafft.fasta"
    output:
        trimmed=f"{PHYLOGENY_DIR}/{{sample}}/Alignment/{{sample}}_SSU_mafft_trimmed.fasta"
    log:
        f"{PHYLOGENY_DIR}/{{sample}}/logs/{{sample}}_trimal.log"
    benchmark:
        f"{PHYLOGENY_DIR}/{{sample}}/benchmarks/{{sample}}_trimal.txt"
    shell:
        r"""
        set -euo pipefail

        trimal \
            -in {input.alignment} \
            -out {output.trimmed} \
            -automated1 \
            > {log} 2>&1
        """
