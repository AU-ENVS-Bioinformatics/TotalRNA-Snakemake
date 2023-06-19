rule infernal_cmsearch:
    input:
        "results/mRNA/trinity.Trinity.fasta"
    output:
        "results/mRNA/trinity/Trinity_cmsearch.out"
    log:
        "logs/cmsearch.log"
    threads: int(config.get("filter_non_coding_rna-THREADS", 50))
    conda:
        "../envs/CoMW.yaml"
    params:
        database= config.get("RFAM_DATABASE"),
    shell:
        "cmsearch -o {output} --cpu {threads} "
        "{params.database} {input} > {log} 2>&1"

rule legacy_tempdir:
    output:
        temp(directory("results/mRNA/trinity/TempFiles"))
    shell: 
        "mkdir -p {output}"

rule parser_cmsearch:
    input:
        tmp = "results/mRNA/trinity/TempFiles",
        infile = "results/mRNA/trinity/Trinity_cmsearch.out"
    output:
        "results/mRNA/trinity/Trinity_cmsearchncRNA.txt"
    log:
        "logs/parser_cmsearch.log"
    conda:
        "../envs/base_python.yaml"
    shell:
        "python workflow/scripts/CoMW/utils/parsecm.py  {input.infile} 1E-3 > {log} 2>&1"

rule exclude_non_coding_rna:
    input:
        fasta = "results/mRNA/trinity.Trinity.fasta",
        exclude = "results/mRNA/trinity/Trinity_cmsearchncRNA.txt"
    output:
        "results/mRNA/contigs_ncrna_filtered.fasta"
    log:
        "logs/exclude_non_coding.log"
    conda:
        "../envs/biopython.yaml"
    script:
        "../scripts/filter_fasta.py"

rule filter_before_align:
    input:
        table="results/mRNA/mapped_reads_to_contigs_AbundanceFiltered.tsv",
        fasta="results/mRNA/trinity/contigs_ncrna_filtered.fasta",
    output:
        "results/mRNA/contigs_ncrna_filtered_AbundanceFiltered.fasta",
    conda:
        "../envs/biopython.yaml"
    log:
        "logs/align_contigs_to_database/filter_abundance.log",
    benchmark:
        "benchmarks/align_contigs_to_database/filter_abundance.log"
    script:
        "../scripts/filter_fasta_by_abundance.py"