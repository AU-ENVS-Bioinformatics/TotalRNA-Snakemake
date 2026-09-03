rule trimal_ssu:
    conda:
        "../envs/trimal.yaml"
    message:
        "[trimAl] Trim the cross-sample SSU alignment"
    input:
        alignment=f"{PHYLOGENY_DIR}/mafft/cross_sample_SSU_reference_mafft.fasta"
    output:
        trimmed=f"{PHYLOGENY_DIR}/mafft/cross_sample_SSU_reference_mafft_trim.fasta"
    log:
        stdout=f"{PHYLOGENY_DIR}/mafft/logs/trimal.log"
    benchmark:
        f"{PHYLOGENY_DIR}/mafft/benchmarks/trimal.txt"
    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {output.trimmed})

        trimal \
            -in {input.alignment} \
            -out {output.trimmed} \
            -automated1 \
            > "{log.stdout}" 2>&1
        """