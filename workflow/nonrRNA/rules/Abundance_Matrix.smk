SALMON_METRICS = [
    "TPM",
    "NumReads",
]

def salmon_quants(wc):
    return expand(
        f"{RESULTS_DIR}/nonrRNA/salmon/{wc.reference}/{{sample}}/quant.sf",
        sample=SAMPLES,
    )

rule salmon_matrix:
    conda:
        "../envs/r_env.yaml"
    message:
        "[Salmon] Build {wildcards.metric} matrix for {wildcards.reference}"
    input:
        quants=salmon_quants
    output:
        abundance_matrix=f"{RESULTS_DIR}/nonrRNA/salmon/{{reference}}/{{metric}}.tsv"
    log:
        stdout=f"{RESULTS_DIR}/nonrRNA/salmon/logs/salmon_matrix_{reference}_{metric}.log"
    benchmark:
        f"{RESULTS_DIR}/nonrRNA/salmon/benchmarks/salmon_matrix_{reference}_{metric}.txt"
    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {output.abundance_matrix})

        Rscript workflow/scripts/build_salmon_matrix.R \
            {wildcards.metric} \
            {output.abundance_matrix} \
            {input.quants} \
            > {log.stdout} 2>&1
        """