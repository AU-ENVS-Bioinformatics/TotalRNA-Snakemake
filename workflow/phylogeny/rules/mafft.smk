################################################################################
# 7. Multiple sequence alignment with MAFFT
################################################################################

rule mafft_ssu:
    conda:
        "../envs/mafft.yaml"
    wildcard_constraints:
        kingdom=BARRNAP_KINGDOM_PATTERN
    message:
        "[MAFFT] Align reconstructed SSUs and SILVA references"
    input:
        fasta=f"{PHYLOGENY_DIR}/{{kingdom}}/references/samples_SSU_{{kingdom}}_with_references.fasta"
    output:
        alignment=f"{PHYLOGENY_DIR}/{{kingdom}}/mafft/samples_SSU_{{kingdom}}_references_mafft.fasta"
    log:
        stdout=f"{PHYLOGENY_DIR}/{{kingdom}}/mafft/logs/mafft.log",
        time=f"{PHYLOGENY_DIR}/{{kingdom}}/mafft/benchmarks/mafft_time.txt"
    benchmark:
        f"{PHYLOGENY_DIR}/{{kingdom}}/mafft/benchmarks/mafft.txt"
    threads:
        config["phylogeny"]["mafft"].get("threads", 8)
    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {output.alignment})

        /usr/bin/time -v -o {log.time} \
            mafft \
                --auto \
                --thread {threads} \
                "{input.fasta}" \
                > "{output.alignment}" \
                2> "{log.stdout}"
        """