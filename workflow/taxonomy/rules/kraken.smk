rule kraken2:
    conda:
        "../envs/kraken2.yaml"
    message:
        "[Kraken2] estimate species composition for {wildcards.sample} for taxonomic classification against {wildcards.database}"
    input:
        rRNA_ssu_r1=f"{RNA_CLASSIFIED_DIR}/{{sample}}/SSU/{{sample}}_SSU_1.fastq.gz",
        rRNA_ssu_r2=f"{RNA_CLASSIFIED_DIR}/{{sample}}/SSU/{{sample}}_SSU_2.fastq.gz",
    output:
        report=f"{TAXONOMY_DIR}/{{sample}}/SSU/{{database}}/kraken/{{sample}}.{{database}}.k2report",
        kraken=f"{TAXONOMY_DIR}/{{sample}}/SSU/{{database}}/kraken/{{sample}}.{{database}}.kraken"
    log:
        stdout=f"{TAXONOMY_DIR}/{{sample}}/SSU/{{database}}/logs/{{sample}}.{{database}}.kraken2.log"
    benchmark:
        f"{TAXONOMY_DIR}/{{sample}}/SSU/{{database}}/benchmarks/{{sample}}.{{database}}.kraken2.txt"
    params:
        db=kraken_db,
        options=config["taxonomy"]["SSU"]["kraken2"].get("options","")
    threads:
        config["taxonomy"]["SSU"]["kraken2"].get("threads", 2)
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.report})

        kraken2 \
          --db {params.db} \
          --threads {threads} \
          --report {output.report} \
          --output {output.kraken} \
          --paired {input.rRNA_ssu_r1} {input.rRNA_ssu_r2} \
            {params.options} \
          > {log.stdout} 2>&1
        """

#kraken2 --db /data_2/Databases/silva_kraken_db/SILVA_138_2_k2db/ --paired AU_RS025_fwd.fq.gz AU_RS025_rev.fq.gz --threads 8 --report kraken2/AU_RS025.report.txt --output kraken2/AU_RS025.output.txt > logs/AU_RS025.kraken2.log