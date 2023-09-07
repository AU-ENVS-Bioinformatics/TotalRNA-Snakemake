trim_galore_params = config.get("trim_galore", "")
trim_galore_threads = config.get("TRIM_GALORE-THREADS", 7)
trim_galore_params.append(f"--cores {trim_galore_threads}")


rule trim_galore_pe:
    input:
        [
            "results/renamed_raw_reads/{sample}_R1.fastq.gz",
            "results/renamed_raw_reads/{sample}_R2.fastq.gz",
        ],
    output:
        fasta_fwd="results/trimmed/{sample}_R1.fq.gz",
        report_fwd="results/trimmed/reports/{sample}_R1_trimming_report.txt",
        fasta_rev="results/trimmed/{sample}_R2.fq.gz",
        report_rev="results/trimmed/reports/{sample}_R2_trimming_report.txt",
    threads: trim_galore_threads
    params:
        extra=" ".join(trim_galore_params),
    log:
        "logs/trim_galore/{sample}.log",
    wrapper:
        "v2.6.0/bio/trim_galore/pe"
