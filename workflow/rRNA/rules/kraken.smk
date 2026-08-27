rule kraken2:
    conda:
        "../envs/kraken2.yaml"
    message:
        "[Kraken2] estimate species composition for {wildcards.sample} for taxonomic classification against {wildcards.database}"
    input:
        rRNA_R1=f"{RESULTS_DIR}/rRNA/{{sample}}/filtered/{{sample}}_rRNA_1.fastq.gz",
        rRNA_R2=f"{RESULTS_DIR}/rRNA/{{sample}}/filtered/{{sample}}_rRNA_2.fastq.gz",
    output:
        report=f"{RESULTS_DIR}/rRNA/{{sample}}/classification/{{database}}/{{sample}}.{{database}}.k2report",
        kraken=f"{RESULTS_DIR}/rRNA/{{sample}}/classification/{{database}}/{{sample}}.{{database}}.kraken"
    log:
        stdout=f"{RESULTS_DIR}/rRNA/{{sample}}/logs/kraken2_{{database}}.log"
    benchmark:
        f"{RESULTS_DIR}/rRNA/{{sample}}/benchmarks/kraken2_{{database}}.txt"
    params:
        db=kraken_db,
        options=config["rRNA"]["kraken2"].get("options","")
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
            {params.options} \
          > {log.stdout} 2>&1
        """

#kraken2 --db /data_2/Databases/silva_kraken_db/SILVA_138_2_k2db/ --paired AU_RS025_fwd.fq.gz AU_RS025_rev.fq.gz --threads 8 --report kraken2/AU_RS025.report.txt --output kraken2/AU_RS025.output.txt > logs/AU_RS025.kraken2.log