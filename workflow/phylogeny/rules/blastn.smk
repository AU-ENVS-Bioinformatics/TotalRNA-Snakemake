################################################################################
# BLASTN PHYLOGENETIC CLASSIFICATION
################################################################################

################################################################################
# BLAST RECONSTRUCTED SSU SEQUENCES AGAINST SILVA
################################################################################

rule blast_cross_sample_ssu:
    conda:
        "../envs/blastn.yaml"
    message:
        "[BLASTN] Identify SILVA references for combined reconstructed SSUs"
    input:
        query=f"{ASSEMBLY_DIR}/concatenated/SSU/cross_sample_SSU.fasta"
    output:
        tsv=f"{PHYLOGENY_DIR}/cross_sample/blastn/cross_sample_SSU_blastn.tsv"
    log:
        stdout=f"{PHYLOGENY_DIR}/cross_sample/logs/cross_sample_blastn.log",
        time=f"{PHYLOGENY_DIR}/cross_sample/benchmarks/cross_sample_blastn_time.txt"
    benchmark:
        f"{PHYLOGENY_DIR}/cross_sample/benchmarks/cross_sample_blastn.txt"
    params:
        db=config["phylogeny"]["blastn"]["database"],
        options=config["phylogeny"]["blastn"]["options"]
    threads:
        config["phylogeny"]["blastn"].get("threads", 12)

    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {output.tsv})

        /usr/bin/time -v -o "{log.time}" \
            blastn \
                -query "{input.query}" \
                -db "{params.db}" \
                -out "{output.tsv}" \
                -num_threads {threads} \
                {params.options} \
                > "{log.stdout}" 2>&1
        """