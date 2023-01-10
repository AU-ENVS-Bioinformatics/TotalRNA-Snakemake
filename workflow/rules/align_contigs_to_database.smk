rule filter_before_align:
    input:
        table="results/mRNA/mapped_reads_to_contigs_AbundanceFiltered.tsv",
        fasta="results/mRNA/trinity/contigs_ncrna_filtered.fasta",
    output:
        "results/mRNA/trinity/contigs_ncrna_filtered_AbundanceFiltered.fasta",
    conda:
        "../envs/biopython.yaml"
    log:
        "logs/align_contigs_to_database/filter_abundance.log",
    benchmark:
        "benchmarks/align_contigs_to_database/filter_abundance.log"
    script:
        "../scripts/filter_fasta_by_abundance.py"


rule translate_trinity:
    input:
        ancient(choose_filepath(config)),
    output:
        f"results/mRNA/trinity/translated/contigs_ncrna_filtered.fasta",
    params:
        orfs=int(config.get("ORFs_translate", 6)),
    conda:
        "../envs/align_contigs_to_database.yaml"
    log:
        "logs/align_contigs_to_database/translate.log",
    benchmark:
        "benchmarks/align_contigs_to_database/translate.log"
    shell:
        "transeq -sformat pearson -clean -frame "
        "{params.orfs} -sequence {input} -outseq {output} "
        ">> {log} 2>&1"


rule split_fasta:
    input:
        f"results/mRNA/trinity/translated/contigs_ncrna_filtered.fasta",
    output:
        expand(
            f"results/mRNA/trinity/translated/contigs_ncrna_filtered.{{index}}.fasta",
            index=parse_pyfasta_int(int(config.get("split_fasta", 2))),
        ),
    params:
        n=int(config.get("split_fasta", 2)),
    conda:
        "../envs/align_contigs_to_database.yaml"
    log:
        "logs/align_contigs_to_database/split.log",
    benchmark:
        "benchmarks/align_contigs_to_database/split.log"
    shell:
        "pyfasta split -n {params.n} {input} "
        ">> {log} 2>&1"


rule sword:
    input:
        f"results/mRNA/trinity/translated/contigs_ncrna_filtered.{{index}}.fasta",
    output:
        f"results/mRNA/trinity/sword/temp/{{database}}.{{index}}.bm9",
    conda:
        "../envs/align_contigs_to_database.yaml"
    params:
        database_path=lambda wildcards: sword_databases.get(wildcards["database"]),
    log:
        f"logs/align_contigs_to_database/sword_{{database}}_{{index}}.log",
    benchmark:
        f"benchmarks/align_contigs_to_database/sword_{{database}}_{{index}}.log"
    threads: int(config.get("align_contigs_to_database-THREADS", 50))
    shell:
        "sword -i {input} "
        "-t {threads} -o {output} "
        "-f bm9 -j {params.database_path} -c 30000 "
        ">> {log} 2>&1"


rule align_contigs_to_database:
    input:
        expand(
            "results/mRNA/trinity/sword/temp/{{database}}.{index}.bm9",
            index=parse_pyfasta_int(int(config.get("split_fasta", 2))),
        ),
    output:
        f"results/mRNA/sword/SWORD_{{database}}_result.bm9",
    conda:
        "../envs/align_contigs_to_database.yaml"
    log:
        f"logs/align_contigs_to_database/cat_{{database}}.log",
    benchmark:
        "benchmarks/align_contigs_to_database/cat_{{database}}.log"
    shell:
        "cat {input} > {output} && "
        "echo Finished > {log}"


rule parse_sword:
    input:
        f"results/mRNA/sword/SWORD_{{database}}_result.bm9",
    output:
        f"results/sword/parsed_{{database}}.tsv",
    conda:
        "../envs/dask.yaml"
    log:
        "logs/sword/parse_sword_{database}.log",
    benchmark:
        "benchmarks/sword/parse_sword_{database}.log"
    script:
        "../scripts/parse_sword.py"
