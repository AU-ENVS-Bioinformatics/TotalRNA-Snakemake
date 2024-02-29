rule count_reads:
    input:
        raw_reads=[
            [
                f"results/renamed_raw_reads/{sample}_R1.fastq.gz"
                for sample in unique_samples
            ],
            [
                f"results/renamed_raw_reads/{sample}_R2.fastq.gz"
                for sample in unique_samples
            ],
        ],
        trimmed_reads=[
            [f"results/trimmed/{sample}_R1.fq.gz" for sample in unique_samples],
            [f"results/trimmed/{sample}_R2.fq.gz" for sample in unique_samples],
        ],
        SSU_sortmerna=[
            [f"results/rrna/{sample}_fwd.fq.gz" for sample in unique_samples],
            [f"results/rrna/{sample}_rev.fq.gz" for sample in unique_samples],
        ],
        not_LSU_sortmerna=[
            [
                f"results/sortmerna/not_LSU/{sample}_fwd.fq.gz"
                for sample in unique_samples
            ],
            [
                f"results/sortmerna/not_LSU/{sample}_rev.fq.gz"
                for sample in unique_samples
            ],
        ],
        rRNA="results/MetaRib/all.dedup.filtered.fasta",
        mRNA="results/mRNA/trinity.Trinity.fasta",
        filtered_mRNA="results/mRNA/abundance_filtered/Trinity_contigs_ncrna_filtered.fasta",
    output:
        "qc/counts/nsequences_file.csv",
    threads: 8
    script:
        "../scripts/qc_read_counts.py"


rule plot_n_sequences:
    input:
        "qc/counts/nsequences_file.csv",
    conda:
        "../envs/phyloseq.yaml"
    output:
        plot1="qc/counts/reads.pdf",
        plot2="qc/counts/sequences_main_steps.pdf",
    script:
        "../scripts/qc_read_counts.R"
