################################################################################
# 7. Multiple sequence alignment with MAFFT
################################################################################

rule mafft_ssu:
    conda:
        "../envs/mafft.yaml"
    message:
        "[MAFFT] align reconstructed SSUs and SILVA references for {wildcards.sample}"
    input:
        fasta=f"{PHYLOGENY_DIR}/{{sample}}/phyloflashreference/{{sample}}_SSU_with_references.fasta"
    output:
        alignment=f"{PHYLOGENY_DIR}/{{sample}}/mafft/{{sample}}_SSU_mafft.fasta"
    log:
        stdout=f"{PHYLOGENY_DIR}/{{sample}}/logs/{{sample}}_mafft.log",
        benchmark_file=f"{PHYLOGENY_DIR}/{{sample}}/benchmarks/{{sample}}_mafft_time.txt",
    benchmark:
        f"{PHYLOGENY_DIR}/{{sample}}/benchmarks/{{sample}}_mafft.txt"
    threads:
        config["phylogeny"]["mafft"].get("threads", 8)
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.alignment})

        /usr/bin/time -v -o {log.benchmark_file} \
            mafft \
                --auto \
                --thread {threads} \
                {input.fasta} \
                > {output.alignment} \
                2> {log.stdout}
        """
