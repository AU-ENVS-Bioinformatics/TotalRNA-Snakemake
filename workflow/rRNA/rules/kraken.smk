rule kraken2:
    conda:
        "../envs/kraken2.yaml"
    message:
        "[Kraken2] estimate species composition for {wildcards.sample} for taxonomic classification"
    input:
        rRNA_R1=f"{RESULTS_DIR}/rRNA/{{sample}}/filtered/{{sample}}_SSU.rRNA.R1.fastq.gz",
        rRNA_R2=f"{RESULTS_DIR}/rRNA/{{sample}}/filtered/{{sample}}_SSU.rRNA.R2.fastq.gz",
    output:
        report=f"{RESULTS_DIR}/rRNA/{{sample}}/classification/{{sample}}.k2report",
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
          --paired {input.rRNA_R1} {input.rRNA_R2} \
          > {log.stdout} 2>&1
        """

#kraken2 --db /data_2/Databases/silva_kraken_db/SILVA_138_2_k2db/ --paired AU_RS025_fwd.fq.gz AU_RS025_rev.fq.gz --threads 8 --report kraken2/AU_RS025.report.txt --output kraken2/AU_RS025.output.txt > logs/AU_RS025.kraken2.log