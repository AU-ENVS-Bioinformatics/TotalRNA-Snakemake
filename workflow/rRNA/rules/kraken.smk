rule kraken2:
    conda:
        "../envs/kraken2.yaml"
    message:
        "[Kraken2] estimate species composition for {wildcards.sample} for taxonomic classification"
    input:
        trim_r1=f"{RESULTS_DIR}/qc/{{sample}}/trimmed/{{sample}}_R1.fastq.gz",
        trim_r2=f"{RESULTS_DIR}/qc/{{sample}}/trimmed/{{sample}}_R2.fastq.gz",
    output:
        report=f"{RESULTS_DIR}/rRNA/{{sample}}/classification/{{sample}}.kreport",
        kraken=f"{RESULTS_DIR}/rRNA/{{sample}}/classification/{{sample}}.kraken"
    log:
        stdout=f"{RESULTS_DIR}/rRNA/{{sample}}/logs/kraken2.log"
    benchmark:
        f"{RESULTS_DIR}/rRNA/{{sample}}/benchmarks/kraken2.txt"
    params:
        db=config["databases"]["kraken_rRNA_db"],
    threads:
        config["rRNA"]["kraken2"].get("threads", 2)
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.report})

        kraken2 \
          --db {params.db} \
          --threads {threads} \
          --report {output.report} \
          --output {output.kraken} \
          --paired {input.trim_r1} {input.trim_r2} \
          > {log.stdout} 2>&1
        """
