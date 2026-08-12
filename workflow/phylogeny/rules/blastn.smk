################################################################################
# 2. BLAST reconstructed SSU sequences against SILVA
################################################################################
    
rule blast_ssu:
    conda:
        "../envs/blastn.yaml"
    message:
        "[BLASTN] identify closest SILVA references for reconstructed SSUs of {wildcards.sample}"
    input:
        phyloflash_query=f"{PHYLOGENY_DIR}/{{sample}}/phyloflash/{{sample}}.all.final.fasta",
    output:
        tsv=f"{PHYLOGENY_DIR}/{{sample}}/blastn/{{sample}}_SSU_blastn.tsv"
    log:
        stdout=f"{PHYLOGENY_DIR}/{{sample}}/logs/{{sample}}_blastn.log",
        benchmark_file=f"{PHYLOGENY_DIR}/{{sample}}/benchmarks/{{sample}}_blastn_time.txt",
    benchmark:
        f"{PHYLOGENY_DIR}/{{sample}}/benchmarks/{{sample}}_blastn.txt"
    params:
        db=config["databases"]["SILVA_138.2_N99_blast"],
        options=config["phylogeny"]["blastn"]["options"]
    threads:
        config["phylogeny"]["blastn"].get("threads", 12)
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.tsv})

        /usr/bin/time -v -o {log.benchmark_file} \
            blastn \
                -query {input.phyloflash_query} \
                -db {params.db} \
                -out {output.tsv} \
                -num_threads {threads} \
                {params.options} \
                > {log.stdout} 2>&1
        """
