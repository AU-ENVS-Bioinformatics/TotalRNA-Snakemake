################################################################################
# BLAST SAMPLE-DERIVED SSUs AGAINST DOMAIN-SPECIFIC REFERENCES
################################################################################

rule blast_samples_SSU_domains:
    conda:
        "../envs/blastn.yaml"
    wildcard_constraints:
        domain=PHYLOGENY_DOMAIN_PATTERN
    message:
        "[BLASTN] Identify SSU references for domain: {wildcards.domain}"
    input:
        query=f"{ASSEMBLY_DIR}/cross_assembly/SSU/domains/{{domain}}/final/transcript.fasta"
    output:
        tsv=f"{PHYLOGENY_DIR}/{{domain}}/blastn/SSU_blastn_{{domain}}.tsv"
    log:
        stdout=f"{PHYLOGENY_DIR}/{{domain}}/blastn/logs/blastn.log",
        time=f"{PHYLOGENY_DIR}/{{domain}}/blastn/benchmarks/blastn_time.txt"
    benchmark:
        f"{PHYLOGENY_DIR}/{{domain}}/blastn/benchmarks/blastn_benchmark.tsv"
    params:
        db=lambda wc: config["phylogeny"]["blastn"]["databases"][wc.domain],
        options=config["phylogeny"]["blastn"].get("options", "")
    threads:
        config["phylogeny"]["blastn"].get("threads", 12)
    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {output.tsv}) $(dirname {log.stdout}) $(dirname {log.time})

        if [ ! -s {input.query} ]; then
            : > {output.tsv}
            : > {log.time}
            printf 'No SSU candidates available for {wildcards.domain}.\n' > {log.stdout}
        else
            /usr/bin/time -v -o {log.time} \
                blastn \
                    -query {input.query} \
                    -db {params.db} \
                    -out {output.tsv} \
                    -num_threads {threads} \
                    {params.options} \
                    > {log.stdout} 2>&1
        fi
        """