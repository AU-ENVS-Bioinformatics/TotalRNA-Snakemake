rule count_reads:
    input:
        raw_reads=[
            [samples_dict[sample].r1 for sample in unique_samples],
            [samples_dict[sample].r2 for sample in unique_samples],
        ],
        trimmed_reads=[
            [f"results/trim_galore/{sample}_R1.fq.gz" for sample in unique_samples],
            [f"results/trim_galore/{sample}_R2.fq.gz" for sample in unique_samples],
        ],
        SSU_sortmerna=[
            [f"results/sortmerna/SSU/{sample}_fwd.fq.gz" for sample in unique_samples],
            [f"results/sortmerna/SSU/{sample}_rev.fq.gz" for sample in unique_samples],
        ],
        not_LSU_sortmerna=[
            [f"results/sortmerna/not_LSU/{sample}_fwd.fq.gz" for sample in unique_samples],
            [f"results/sortmerna/not_LSU/{sample}_rev.fq.gz" for sample in unique_samples],
        ],
        rRNA="results/metarib/final_contigs.fasta",
        mRNA="results/trinity/trinity.Trinity.fasta",
        filtered_mRNA="results/mRNA/Trinity_contigs_cmsearch_filtered.fasta",
    conda:
        "../envs/base_python.yaml"
    output:
        "qc/counts/nsequences_file.csv",
    log:
        "logs/qc/count_sequences.log",
    threads: 8
    script:
        "../scripts/qc_read_counts.py"


rule plot_n_sequences:
    input:
        "qc/counts/nsequences_file.csv",
    conda:
        "../envs/phyloseq.yaml"
    log:
        "logs/qc/plot_count_sequences.log",
    output:
        plot1="qc/counts/reads.pdf",
        plot2="qc/counts/sequences_main_steps.pdf",
    script:
        "../scripts/qc_read_counts.R"
