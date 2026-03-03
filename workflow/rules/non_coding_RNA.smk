rule infernal_cmsearch:
    input:
        fasta="results/trinity/trinity.Trinity.fasta",
        database=config.get("RFAM_DATABASE"),
    output:
        out="results/cmsearch/RFAM_cmsearch.out",
        tbl="results/cmsearch/RFAM_cmsearch.tbl",
    log:
        "logs/cmsearch.log",
    threads: config["threads"]["cmsearch"]
    conda:
        "../envs/infernal.yaml"
    shell:
        "cmsearch --tblout {output.tbl} -o {output.out} --cpu {threads} "
        "{input.database} {input.fasta} > {log} 2>&1"


rule processing_cmsearch_tbl:
    input:
        "results/cmsearch/RFAM_cmsearch.tbl",
    output:
        "results/cmsearch/RFAM_cmsearch.tbl.processed.tsv",
    log:
        "logs/processing_cmsearch_tbl.log",
    conda:
        "../envs/pandas.yaml"
    params:
        evalue=config["evalue_non_coding_rna"],
    script:
        "../scripts/processing_cmsearch_tbl.py"


rule exclude_non_coding_rna:
    input:
        fasta="results/trinity/trinity.Trinity.fasta",
        exclude="results/cmsearch/RFAM_cmsearch.tbl.processed.tsv",
    output:
        "results/mRNA/Trinity_contigs_cmsearch_filtered.fasta",
    log:
        "logs/exclude_non_coding.log",
    conda:
        "../envs/biopython.yaml"
    script:
        "../scripts/filter_fasta.py"