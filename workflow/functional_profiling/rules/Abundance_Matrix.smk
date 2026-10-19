rule salmon_matrix:
    conda:
        "../envs/r_env.yaml"
    wildcard_constraints:
        reference=SALMON_REFERENCE_PATTERN,
        metric=SALMON_METRIC_PATTERN
    message:
        "[Salmon] Build {wildcards.metric} matrix for {wildcards.reference}"
    input:
        quants=lambda wc: expand(
            f"{FUNCTION_DIR}/salmon/{wc.reference}/{{sample}}/quant.sf",
            sample=SAMPLES,
        )
    output:
        matrix=f"{FUNCTION_DIR}/salmon/{{reference}}/{{metric}}.tsv"
    log:
        stdout=f"{FUNCTION_DIR}/salmon/logs/salmon_matrix_{{reference}}_{{metric}}.log"
    benchmark:
        f"{FUNCTION_DIR}/salmon/benchmarks/salmon_matrix_{{reference}}_{{metric}}.txt"
    params:
        script=f"{FUNCTIONAL_SCRIPTS_DIR}/build_salmon_matrix.R"
    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {output.matrix})

        Rscript {params.script} \
            {wildcards.metric} \
            {output.matrix} \
            {input.quants} \
            > {log.stdout} 2>&1
        """