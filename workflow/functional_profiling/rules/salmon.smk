################################################################################
# SALMON QUANTIFICATION
################################################################################

FUNCTIONAL_MODULE = config["functional_profiling"]["module"]

SALMON_REFERENCE_FASTAS = {
    "transcripts": f"{ASSEMBLY_DIR}/coassembly/non_rRNA/rnaspades/transcripts.fasta",
    "cds": f"{FUNCTION_DIR}/ORF_prediction/final_orf/ORF_genes.fasta",
}

SALMON_REFERENCE_PATTERN = "|".join(SALMON_REFERENCE_FASTAS)


def salmon_reference(wildcards):
    return SALMON_REFERENCE_FASTAS[wildcards.reference]


if FUNCTIONAL_MODULE in ["assembly", "both", "all"]:

    ################################################################################
    # BUILD SALMON INDEX
    ################################################################################

    rule salmon_index:
        conda:
            "../envs/salmon.yaml"
        wildcard_constraints:
            reference=SALMON_REFERENCE_PATTERN
        message:
            "[Salmon] Build index for {wildcards.reference}"
        input:
            reference=salmon_reference
        output:
            index=directory(f"{FUNCTION_DIR}/salmon/{{reference}}_index")
        log:
            stdout=f"{FUNCTION_DIR}/salmon/logs/salmon_index_{{reference}}.log"
        benchmark:
            f"{FUNCTION_DIR}/salmon/benchmarks/salmon_index_{{reference}}.txt"
        threads:
            config["functional_profiling"]["salmon"].get("threads", 4)
        shell:
            r"""
            set -euo pipefail

            mkdir -p $(dirname {output.index})

            salmon index \
                -t {input.reference} \
                -i {output.index} \
                -p {threads} \
                > {log.stdout} 2>&1
            """

    ################################################################################
    # QUANTIFY EACH SAMPLE
    ################################################################################

    rule salmon_quant:
        conda:
            "../envs/salmon.yaml"
        wildcard_constraints:
            reference=SALMON_REFERENCE_PATTERN
        message:
            "[Salmon] Quantify {wildcards.sample} against {wildcards.reference}"
        input:
            index=f"{FUNCTION_DIR}/salmon/{{reference}}_index",
            r1=f"{RNA_CLASSIFIED_DIR}/{{sample}}/non_rRNA/{{sample}}_non_rRNA_1.fastq.gz",
            r2=f"{RNA_CLASSIFIED_DIR}/{{sample}}/non_rRNA/{{sample}}_non_rRNA_2.fastq.gz"
        output:
            quant=f"{FUNCTION_DIR}/salmon/{{reference}}/{{sample}}/quant.sf"
        log:
            stdout=f"{FUNCTION_DIR}/salmon/logs/{{reference}}_{{sample}}.log"
        benchmark:
            f"{FUNCTION_DIR}/salmon/benchmarks/{{reference}}_{{sample}}.txt"
        params:
            options=config["functional_profiling"]["salmon"].get("options", "")
        threads:
            config["functional_profiling"]["salmon"].get("threads", 4)
        shell:
            r"""
            set -euo pipefail

            mkdir -p $(dirname {output.quant})
            
            salmon quant \
                -i {input.index} \
                -l A \
                -1 {input.r1} \
                -2 {input.r2} \
                -p {threads} \
                -o $(dirname {output.quant}) \
                {params.options} \
                > {log.stdout} 2>&1
            """