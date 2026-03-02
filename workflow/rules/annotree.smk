rule diamond:
    input:
        fasta="results/trinity/trinity.Trinity.fasta",
        database=config["ANNOTREE"]["DATABASE"],
    output:
        "results/diamond/annotree.tsv",
    log:
        "logs/diamond/run_diamond.log",
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


rule export_mrna_abundance_tsv:
    input:
        db="results/mRNA/database.db",
    output:
        "results/mRNA/mapped_reads_to_contigs.tsv",
    params:
        table="abundance",
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/mRNA/export_abundance_mrna.log",
    script:
        "../scripts/abundance2tsv.py"


rule annotate_and_func_ann:
    input:
        tsv="results/diamond/annotree.tsv",
        mapping=config["ANNOTREE"]["MAPPING"],
        counts="results/mRNA/mapped_reads_to_contigs.tsv",
        brite=config["BRITE"],
    output:
        annotated="results/diamond/annotree_ann.tsv",
        brite="results/diamond/annotree_brite.tsv",
        meta="results/diamond/annotree_meta.tsv",
    log:
        "logs/diamond/annotate_and_func_ann.log",
    params:
        min_score_threshold=100,
        threshold=0.95,
    conda:
        "../envs/annotree.yaml"
    script:
        "../scripts/annotate_and_func_ann.py"


# filtered="results/mRNA/filter_contigs.done",
# "awk '{{print $1}}' {input.filtered} | seqtk subseq {input.fasta} - > temp.fasta 2> {log}"
