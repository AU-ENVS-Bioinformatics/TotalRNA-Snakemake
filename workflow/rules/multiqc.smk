def multiqc_input(wildcards):
    branch = {
        "trim_galore": expand(
            "results/trim_galore/reports/{sample}_R{dir}_trimming_report.txt",
            sample=unique_samples,
            dir=["1", "2"],
        ),
        "sortmerna_SSU": expand(
            "results/sortmerna/SSU/{sample}.log",
            sample=unique_samples,
        ),
        "sortmerna_LSU": expand(
            "results/sortmerna/LSU/{sample}.log",
            sample=unique_samples,
        ),
    }
    return branch[wildcards.qc_type]


rule multiqc:
    input:
        multiqc_input,
    output:
        "qc/{qc_type}_multiqc.html",
    priority: 50
    params:
        extra="",
    log:
        "logs/multiqc/{qc_type}_multiqc.log",
    wrapper:
        "v2.7.0/bio/multiqc"
