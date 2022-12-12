rule translate_filtered_trinity:
    input:
        ancient(f"results/mRNA/trinity/contigs_ncrna_filtered.fasta"),
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
            index=list(range(int(config.get("split_fasta", 2)))),
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
        f"results/mRNA/trinity/sword/temp/{{database}}.{{index}}.tsv",
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
            "results/mRNA/trinity/sword/temp/{{database}}.{index}.tsv",
            index=list(range(int(config.get("split_fasta", 2)))),
        ),
    output:
        f"results/mRNA/ML_SWORD_{{database}}_result.tsv",
    conda:
        "../envs/align_contigs_to_database.yaml"
    log:
        f"logs/align_contigs_to_database/cat_{{database}}.log",
    benchmark:
        "benchmarks/align_contigs_to_database/cat_{{database}}.log"
    shell:
        "cat {input} > {output} && "
        "echo Finished > {log}"
