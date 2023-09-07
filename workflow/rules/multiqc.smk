rule multiqc_dir:
    input:
        expand("results/trimmed/reports/{sample}_R{dir}_trimming_report.txt", sample=unique_samples, dir = ["1", "2"])
    output:
        "qc/trim_galore_multiqc.html"
    params:
        extra=""  # Optional: extra parameters for multiqc.
    log:
        "logs/multiqc.log"
    wrapper:
        "v2.6.0/bio/multiqc"