################################################################################
# 8. Trim alignment with trimAl
################################################################################

rule trimal_ssu:
    conda:
        "../envs/trimal.yaml"
    message:
        "[trimAl] trim SSU multiple-sequence-alignment (MSAs) across samples and references"
    input:
        alignment=f"{PHYLOGENY_DIR}/mafft/cross_sample_SSU_reference_mafft.fasta"
    output:
        trimmed=f"{PHYLOGENY_DIR}/mafft/cross_sample_SSU_reference_mafft_trim.fasta"
    log:
        f"{PHYLOGENY_DIR}/mafft/logs/trimal.log"
    benchmark:
        f"{PHYLOGENY_DIR}/mafft/benchmarks/trimal.txt"
    shell:
        r"""
        set -euo pipefail

        trimal \
            -in {input.alignment} \
            -out {output.trimmed} \
            -automated1 \
            > {log} 2>&1
        """
