rule diamond:
    input:
        fasta="results/trinity/trinity.Trinity.fasta",
        database=config["ANNOTREE"]["DATABASE"],
    output:
        "results/annotree/annotree.tsv",
    log:
        "logs/annotree/run_diamond.log",
    threads: config["threads"]["diamond"]
    conda:
        "../envs/annotree.yaml"
    shell:
        "printf 'qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\n' > {output} && "
        "diamond blastx "
        "--query {input.fasta} "
        "--db {input.database} "
        "--threads {threads} "
        "--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore "
        ">> {output} 2> {log}"


rule annotree:
    input:
        tsv="results/annotree/annotree.tsv",
        mapping=config["ANNOTREE"]["MAPPING"],
        counts="results/mRNA/mapped_reads_to_contigs.tsv",
        brite=config["BRITE"],
    output:
        annotated="results/annotree/annotated.tsv",
        brite="results/annotree/brite.tsv",
        meta="results/annotree/metabolism.tsv",
    log:
        "logs/annotree/annotate_and_func_ann.log",
    params:
        min_score_threshold=100,
        threshold=0.95,
    conda:
        "../envs/annotree.yaml"
    script:
        "../scripts/annotate_and_func_ann.py"
