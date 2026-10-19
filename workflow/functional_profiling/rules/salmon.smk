rule salmon_index:
    conda:
        "../envs/salmon.yaml"
    wildcard_constraints:
        reference=SALMON_REFERENCE_PATTERN
    message:
        "[Salmon] Build index for {wildcards.reference}"
    input:
        reference=lambda wc: salmon_reference(
            wc,
            results_dir=RESULTS_DIR,
        )
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