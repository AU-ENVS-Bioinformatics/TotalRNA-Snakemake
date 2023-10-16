rule multiqc_trim_galore:
    input:
        expand(
            "results/trimmed/reports/{sample}_R{dir}_trimming_report.txt",
            sample=unique_samples,
            dir=["1", "2"],
        ),
    output:
        "qc/trim_galore_multiqc.html",
    priority: 50
    params:
        extra="",  # Optional: extra parameters for multiqc.
    log:
        "logs/multiqc/multiqc_trim_galore.log",
    wrapper:
        "v2.6.0/bio/multiqc"


rule multiqc_sortmerna_SSU:
    input:
        expand(
            "results/rrna/{sample}.log",
            sample=unique_samples,
        ),
    output:
        "qc/sortmerna_SSU_multiqc.html",
    priority: 50
    params:
        extra="",  # Optional: extra parameters for multiqc.
    log:
        "logs/multiqc/sortmerna_SSU_multiqc.log",
    wrapper:
        "v2.6.0/bio/multiqc"


rule multiqc_sortmerna_LSU:
    input:
        expand(
            "results/sortmerna/LSU/{sample}.log",
            sample=unique_samples,
        ),
    output:
        "qc/sortmerna_LSU_multiqc.html",
    priority: 50
    params:
        extra="",  # Optional: extra parameters for multiqc.
    log:
        "logs/multiqc/sortmerna_LSU_multiqc.log",
    wrapper:
        "v2.6.0/bio/multiqc"
