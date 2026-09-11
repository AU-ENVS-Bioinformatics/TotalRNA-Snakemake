################################################################################
# SALMON ABUNDANCE MATRICES
################################################################################

FUNCTIONAL_MODULE = config["functional_profiling"]["module"]

SALMON_MATRIX_REFERENCES = [
    "transcripts",
    "cds",
]

SALMON_MATRIX_METRICS = [
    "TPM",
    "NumReads",
]

SALMON_MATRIX_REFERENCE_PATTERN = "|".join(SALMON_MATRIX_REFERENCES)
SALMON_MATRIX_METRIC_PATTERN = "|".join(SALMON_MATRIX_METRICS)

FUNCTIONAL_SCRIPTS_DIR = f"{workflow.basedir}/nonrRNA/scripts"


def salmon_quants(wildcards):
    return expand(
        f"{FUNCTION_DIR}/salmon/{wildcards.reference}/{{sample}}/quant.sf",
        sample=SAMPLES,
    )


if FUNCTIONAL_MODULE in ["assembly", "both", "all"]:

    rule salmon_matrix:
        conda:
            "../envs/r_env.yaml"
        wildcard_constraints:
            reference=SALMON_MATRIX_REFERENCE_PATTERN,
            metric=SALMON_MATRIX_METRIC_PATTERN
        message:
            "[Salmon] Build {wildcards.metric} matrix for {wildcards.reference}"
        input:
            quants=salmon_quants
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