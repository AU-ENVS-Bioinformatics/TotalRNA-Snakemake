rule infernal_cmsearch:
    input:
        fasta="results/mRNA/trinity.Trinity.fasta",
        database=config.get("RFAM_DATABASE"),
    output:
        "results/mRNA/non_coding_rna/RFAM_cmsearch.out",
    log:
        "logs/cmsearch.log",
    threads: int(config.get("filter_non_coding_rna-THREADS", 50))
    conda:
        "../envs/infernal.yaml"
    shell:
        "cmsearch -o {output} --cpu {threads} "
        "{input.database} {input.fasta} > {log} 2>&1"
