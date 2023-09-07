rule multiqc_trim_galore:
    input:
        expand(
            "results/trimmed/reports/{sample}_R{dir}_trimming_report.txt",
            sample=unique_samples,
            dir=["1", "2"],
        ),
    output:
        "qc/trim_galore_multiqc.html",
    params:
        extra="",  # Optional: extra parameters for multiqc.
    log:
        "logs/multiqc_trim_galore.log",
    wrapper:
        "v2.6.0/bio/multiqc"


rule multiqc_sortmerna:
    input:
        expand(
            "results/rrna/{sample}.log",
            sample=unique_samples,
        ),
    output:
        "qc/sortmerna_multiqc.html",
    params:
        extra="",  # Optional: extra parameters for multiqc.
    log:
        "logs/sortmerna_multiqc.log",
    wrapper:
        "v2.6.0/bio/multiqc"
