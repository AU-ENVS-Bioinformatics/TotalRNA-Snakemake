rule infernal_cmsearch:
    input:
        fasta="results/mRNA/trinity.Trinity.fasta",
        database=config.get("RFAM_DATABASE"),
    output:
        out="results/mRNA/non_coding_rna/RFAM_cmsearch.out",
        tbl="results/mRNA/non_coding_rna/RFAM_cmsearch.tbl",
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
        "results/mRNA/non_coding_rna/RFAM_cmsearch.tbl",
    output:
        "results/mRNA/non_coding_rna/RFAM_cmsearch.tbl.processed.tsv",
    log:
        "logs/processing_cmsearch_tbl.log",
    conda:
        "../envs/pandas.yaml"
    params:
        evalue=config["evalue_non_coding_rna"],
    script:
        "../scripts/processing_cmsearch_tbl.py"


rule non_coding_fasta_names:
    input:
        "results/mRNA/non_coding_rna/RFAM_cmsearch.tbl.processed.tsv",
    output:
        "results/mRNA/non_coding_rna/RFAM_cmsearch_names.txt",
    log:
        "logs/non_coding_headers.log",
    conda:
        "../envs/base_python.yaml"
    shell:
        # Remove first line and print first column
        "grep 'target_name' -v < {input} | cut -f 1 > {output} 2> {log}"


rule exclude_non_coding_rna:
    input:
        fasta="results/mRNA/trinity.Trinity.fasta",
        exclude="results/mRNA/non_coding_rna/RFAM_cmsearch_names.txt",
    output:
        "results/mRNA/Trinity_contigs_ncrna_filtered.fasta",
    log:
        "logs/exclude_non_coding.log",
    conda:
        "../envs/biopython.yaml"
    script:
        "../scripts/filter_fasta.py"
