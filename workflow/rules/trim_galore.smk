trim_galore_params = config.get("trim_galore", "")
trim_galore_threads = min(AVAILABLE_THREADS, 7)
trim_galore_params.append(f"--cores {trim_galore_threads}")


rule trim_galore:
    input:
        [
            ancient(f"results/renamed_raw_reads/{{sample}}_R1.fastq.gz"),
            ancient(f"results/renamed_raw_reads/{{sample}}_R2.fastq.gz"),
        ],
    output:
        f"results/trimmed/{{sample}}_R1_val_1.fq.gz",
        report(
            f"results/trimmed/{{sample}}_R1.fastq.gz_trimming_report.txt",
            caption="report/trim_galore.rst",
            category="Quality Control and Trimming",
            subcategory="{sample}",
        ),
        f"results/trimmed/{{sample}}_R2_val_2.fq.gz",
        report(
            f"results/trimmed/{{sample}}_R2.fastq.gz_trimming_report.txt",
            caption="report/trim_galore.rst",
            category="Quality Control and Trimming",
            subcategory="{sample}",
        ),
    params:
        extra=" ".join(trim_galore_params),
    log:
        "logs/trim_galore/{sample}.log",
    wrapper:
        "v1.12.1/bio/trim_galore/pe"
