trim_galore_params = config.get("trim_galore", "")
trim_galore_threads = config["threads"]["trim_galore"]
trim_galore_params.append(f"--cores {trim_galore_threads}")


rule trim_galore_pe:
    input:
        [
            lambda wc: samples_dict[wc.sample].r1,
            lambda wc: samples_dict[wc.sample].r2,
        ],
    output:
        fasta_fwd="results/trim_galore/{sample}_R1.fq.gz",
        report_fwd="results/trim_galore/reports/{sample}_R1_trimming_report.txt",
        fasta_rev="results/trim_galore/{sample}_R2.fq.gz",
        report_rev="results/trim_galore/reports/{sample}_R2_trimming_report.txt",
    threads: trim_galore_threads
    params:
        extra=" ".join(trim_galore_params),
    log:
        "logs/trim_galore/{sample}.log",
    wrapper:
        "v2.6.0/bio/trim_galore/pe"
