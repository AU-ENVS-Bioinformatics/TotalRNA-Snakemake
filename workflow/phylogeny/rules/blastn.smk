################################################################################
# BLASTN PHYLOGENETIC CLASSIFICATION
################################################################################

################################################################################
# 2. BLAST reconstructed SSU sequences against SILVA across sample
################################################################################
    
rule blast_cross_sample_ssu:
    conda:
        "../envs/blastn.yaml"
    message:
        "[BLASTN] identify SILVA references for combined reconstructed SSUs"
    input:
        query=f"{PHYLOGENY_DIR}/cross_sample/reconstructed_SSU_all_samples.fasta"
    output:
        tsv=f"{PHYLOGENY_DIR}/cross_sample/blastn/cross_sample_SSU_blastn.tsv"
    log:
        stdout=f"{PHYLOGENY_DIR}/cross_sample/logs/cross_sample_blastn.log",
        benchmark_file=f"{PHYLOGENY_DIR}/cross_sample/benchmarks/cross_sample_blastn_time.txt"
    benchmark:
        f"{PHYLOGENY_DIR}/cross_sample/benchmarks/cross_sample_blastn.txt"
    params:
        db=config["databases"]["SILVA_138.2_N99_blast"],
        options=config["phylogeny"]["blastn"]["options"]
    threads:
        config["phylogeny"]["blastn"].get("threads", 12)
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.tsv})

        /usr/bin/time -v -o "{log.benchmark_file}" \
            blastn \
                -query "{input.query}" \
                -db "{params.db}" \
                -out "{output.tsv}" \
                -num_threads {threads} \
                {params.options} \
                > {log.stdout} 2>&1
        """

################################################################################
# 2. BLAST EXTRACTED EUKARYOTIC 18S SEQUENCES AGAINST PR2
################################################################################

rule blast_eukaryotic_ssu_pr2:
    conda:
        "../envs/blastn.yaml"
    message:
        "[BLASTN: PR2] classify reconstructed eukaryotic 18S sequences"
    input:
        query=f"{PHYLOGENY_DIR}/cross_sample/PR2/queries/eukaryota_18S.fasta"
    output:
        tsv=f"{PHYLOGENY_DIR}/cross_sample/PR2/blastn/eukaryota_18S_vs_PR2.tsv"
    log:
        stdout=f"{PHYLOGENY_DIR}/cross_sample/PR2/logs/eukaryota_18S_vs_PR2.log",
        benchmark_file=f"{PHYLOGENY_DIR}/cross_sample/PR2/benchmarks/eukaryota_18S_vs_PR2_time.txt"
    benchmark:
        f"{PHYLOGENY_DIR}/cross_sample/PR2/benchmarks/eukaryota_18S_vs_PR2.txt"
    params:
        db=config["databases"]["PR2_SSU_blast"],
        options=config["phylogeny"]["blastn"]["PR2"]["options"]
    threads:
        config["phylogeny"]["blastn"]["PR2"].get("threads", 12)
    shell:
        r"""
        set -euo pipefail

        mkdir -p \
            "$(dirname "{output.tsv}")" \
            "$(dirname "{log.stdout}")" \
            "$(dirname "{log.benchmark_file}")"

        /usr/bin/time -v -o "{log.benchmark_file}" \
            blastn \
                -query "{input.query}" \
                -db "{params.db}" \
                -out "{output.tsv}" \
                -num_threads {threads} \
                {params.options} \
                > "{log.stdout}" 2>&1
        """