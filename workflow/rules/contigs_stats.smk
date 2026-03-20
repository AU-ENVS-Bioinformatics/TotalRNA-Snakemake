_bwa_index_ext = [".amb", ".ann", ".bwt", ".pac", ".sa"]


rule bwa_index:
    input:
        fasta="results/{folder}/{file}.fasta",
    output:
        idx=multiext("results/{folder}/{file}.fasta", *_bwa_index_ext),
    log:
        "logs/bwa/index/{folder}/{file}.log",
    params:
        algorithm="bwtsw",
    wrapper:
        "v2.6.0/bio/bwa/index"


_contig_reads_config = {
    "rRNA": {
        "contig": "metarib/final_contigs.fasta",
        "read_dir": "sortmerna/SSU",
    },
    "mRNA": {
        "contig": "trinity/trinity.Trinity.fasta",
        "read_dir": "sortmerna/not_LSU",
    },
}


rule bwa_map_and_sort:
    input:
        contig=lambda wc: multiext(
            f"results/{_contig_reads_config[wc.folder]['contig']}", "", *_bwa_index_ext
        ),
        r1=lambda wc: f"results/{_contig_reads_config[wc.folder]['read_dir']}/{wc.sample}_fwd.fq.gz",
        r2=lambda wc: f"results/{_contig_reads_config[wc.folder]['read_dir']}/{wc.sample}_rev.fq.gz",
    output:
        "results/{folder}/bwa/{sample}_sorted.bam",
    log:
        "logs/bwa/{folder}_{sample}.log",
    threads: config["threads"]["bwamem"]
    conda:
        "../envs/bwa_samtools.yaml"
    shell:
        """
        bwa mem -t {threads} {input.contig[0]} {input.r1} {input.r2} 2> {log} \
        | samtools sort -o {output}
        """


rule samtools_idxstats_contigs:
    input:
        "results/{folder}/bwa/{sample}_sorted.bam",
    output:
        idx="results/{folder}/bwa/{sample}_sorted.bam.bai",
        idxstats="results/{folder}/bwa/{sample}_sorted.bam.idxstats",
    log:
        "logs/samtools/idxstats/{folder}_{sample}.log",
    conda:
        "../envs/bwa_samtools.yaml"
    shell:
        """
        samtools index {input} 2> {log}
        samtools idxstats {input} > {output.idxstats} 2>> {log}
        """





rule create_database:
    input:
        idxstats=lambda wc: expand(
            f"results/{wc.folder}/bwa/{{sample}}_sorted.bam.idxstats",
            sample=unique_samples,
        ),
    output:
        database="results/{folder}/database.db",
    log:
        "logs/{folder}/create_database.log",
    params:
        dir_path=lambda wc: f"results/{wc.folder}/bwa",
    conda:
        "../envs/duckdb.yaml"
    script:
        "../scripts/create_database_table.py"


rule mapped_read_length:
    input:
        bam=lambda wc: expand(f"results/{wc.folder}/bwa/{{sample}}_sorted.bam", sample=unique_samples),
        database="results/{folder}/database.db",
    output:
        touch("results/{folder}/read_length.done"),
    log:
        "logs/{folder}/bwa/read_length.log",
    conda:
        "../envs/pysam.yaml"
    script:
        "../scripts/get_read_length.py"


rule export_tsv:
    input:
        db="results/{folder}/database.db",
        reads_done="results/{folder}/read_length.done",
    output:
        "results/{folder}/mapped_reads_exported.tsv",
        "results/{folder}/read_length_exported.tsv",
        "results/{folder}/contig_length_exported.tsv"
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/{folder}/export_tsv.log",
    script:
        "../scripts/export_tsv.py"


# filtered="results/mRNA/filter_contigs.done",
# "awk '{{print $1}}' {input.filtered} | seqtk subseq {input.fasta} - > temp.fasta 2> {log}"
