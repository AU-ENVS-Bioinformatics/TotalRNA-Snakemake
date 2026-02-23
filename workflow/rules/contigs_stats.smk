_bwa_index_ext = ["", ".amb", ".ann", ".bwt", ".pac", ".sa"]


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
        "contig": "metarib/all.dedup.filtered.fasta",
        "read_dir": "sortmerna/SSU",
    },
    "mRNA": {
        "contig": "trinity/trinity.Trinity.fasta",
        "read_dir": "sortmerna/not_LSU",
    },
}


rule bwa_map_to_contigs:
    input:
        contig=lambda wc: multiext(
            f"results/{_contig_reads_config[wc.folder]['contig']}", *_bwa_index_ext
        ),
        read=lambda wc: f"results/{_contig_reads_config[wc.folder]['read_dir']}/{wc.sample_dir}.fq.gz",
    output:
        temp("results/{folder}/bwa/{sample_dir}.sam"),
    log:
        "logs/bwa/{folder}_{sample_dir}.log",
    threads: config["threads"]["bwamem"]
    conda:
        "../envs/bwa.yaml"
    shell:
        """
        bwa mem -t {threads} {input.contig[0]} {input.read} > {output} 2> {log}
        """


rule sort_bwa_contigs:
    input:
        "results/{folder}/bwa/{sample_dir}.sam",
    output:
        "results/{folder}/bwa/{sample_dir}_sorted.bam",
    log:
        "logs/samtools/{folder}_{sample_dir}_sort.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools view -b -S {input} | samtools sort -o {output} - 2> {log}
        """


rule samtools_idxstats_contigs:
    input:
        bam="results/{folder}/bwa/{sample}_{dir}_sorted.bam",
    output:
        idx="results/{folder}/bwa/{sample}_{dir}_sorted.bam.bai",
        idxstats="results/{folder}/bwa/{sample}_{dir}_sorted.bam.idxstats",
        idxstats_tsv="results/{folder}/bwa/{sample}_{dir}_sorted.bam.idxstats.tsv",
    log:
        "logs/samtools/idxstats/{folder}_{sample}_{dir}.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools index {input.bam} 2> {log}
        samtools idxstats {input.bam} > {output.idxstats} 2>> {log}
        awk -v sample="{wildcards.sample}" 'BEGIN{{OFS="\\t"; print "sample", "contig", "mapped_reads"}} {{print sample, $1, $3}}' {output.idxstats} > {output.idxstats_tsv}
        """


rule sample_mapped_read_length:
    input:
        bam="results/{folder}/bwa/{sample_dir}_sorted.bam",
    output:
        tsv="results/{folder}/bwa/{sample_dir}_reads.tsv",
    log:
        "logs/{folder}/bwa/{sample_dir}_reads.log",
    conda:
        "../envs/pysam.yaml"
    script:
        "../scripts/sample_mapped_read_length.py"


rule create_database:
    input:
        single_sample=lambda wc: f"results/{wc.folder}/bwa/{unique_samples[0]}_fwd_sorted.bam.idxstats",
        counts_length_tsv=lambda wc: expand(
            [
                "results/{folder}/bwa/{sample}_{dir}_sorted.bam.idxstats.tsv",
                "results/{folder}/bwa/{sample}_{dir}_reads.tsv",
            ],
            folder=wc.folder,
            sample=unique_samples,
            dir=["rev", "fwd"],
        ),
    output:
        database="results/{folder}/database.db",
    log:
        "logs/{folder}/create_database.log",
    conda:
        "../envs/duckdb.yaml"
    script:
        "../scripts/create_database_table.py"
