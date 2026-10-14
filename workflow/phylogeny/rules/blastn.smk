################################################################################
# BLAST RECONSTRUCTED SSUs AGAINST KINGDOM-SPECIFIC REFERENCES
################################################################################

rule blast_samples_SSU_kingdoms:
    conda:
        "../envs/blastn.yaml"
    wildcard_constraints:
        kingdom=BARRNAP_KINGDOM_PATTERN
    message:
        "[BLASTN] Identify SSU references for kingdom: {wildcards.kingdom}"
    input:
        query=f"{ASSEMBLY_DIR}/SSU_candidates/postprocessing/{{kingdom}}/SSU_{{kingdom}}_catalogue.fasta"
    output:
        tsv=f"{PHYLOGENY_DIR}/{{kingdom}}/blastn/SSU_blastn_{{kingdom}}.tsv"
    log:
        stdout=f"{PHYLOGENY_DIR}/{{kingdom}}/blastn/logs/blastn.log",
        time=f"{PHYLOGENY_DIR}/{{kingdom}}/blastn/benchmarks/blastn_time.txt"
    benchmark:
        f"{PHYLOGENY_DIR}/{{kingdom}}/blastn/benchmarks/blastn_benchmark.tsv"
    params:
        db=lambda wc: config["phylogeny"]["blastn"]["databases"][wc.kingdom],
        options=config["phylogeny"]["blastn"].get("options", "")
    threads:
        config["phylogeny"]["blastn"].get("threads", 12)
    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {output.tsv})

        /usr/bin/time -v -o {log.time} \
            blastn \
                -query {input.query} \
                -db {params.db} \
                -out {output.tsv} \
                -num_threads {threads} \
                {params.options} \
                > {log.stdout} 2>&1
        """