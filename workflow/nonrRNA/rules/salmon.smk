SALMON_REFERENCES = {
    "transcripts": (
        f"{RESULTS_DIR}/nonrRNA/coassembly/transcripts.fasta"
    ),
    "cds": (
        f"{RESULTS_DIR}/nonrRNA/predicted/transcripts.fasta.transdecoder.cds"
    ),
}

def salmon_reference(wc):
    return SALMON_REFERENCES[wc.reference]

rule salmon_index:
    conda:
        "../envs/salmon.yaml"
    message:
        "[Salmon] Build {wildcards.reference} index"
    input:
        transcripts=salmon_reference
    output:
        index=directory(
            f"{RESULTS_DIR}/nonrRNA/salmon/{{reference}}_index"
        )
    log:
        stdout=f"{RESULTS_DIR}/nonrRNA/salmon/logs/salmon_index_{{reference}}.log"
    benchmark:
        f"{RESULTS_DIR}/nonrRNA/salmon/benchmarks/salmon_index_{{reference}}.txt"
    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {output.index})

        salmon index \
            -t {input.transcripts} \
            -i {output.index} \
            > {log.stdout} 2>&1
        """

nonrrna_module = config["nonrRNA"]["nonrrna_module"]

if nonrrna_module in ["coassembly", "both", "all"]:
    rule salmon:
        conda:
            "../envs/salmon.yaml"
        message:
            "[Salmon] Quantifying {wildcards.sample} against {wildcards.reference}"
        input:
            index=f"{RESULTS_DIR}/nonrRNA/salmon/{{reference}}_index",
            nonrRNA_R1=f"{RESULTS_DIR}/nonrRNA/{{sample}}/filtered/{{sample}}_nonrRNA_1.fastq.gz",
            nonrRNA_R2=f"{RESULTS_DIR}/nonrRNA/{{sample}}/filtered/{{sample}}_nonrRNA_2.fastq.gz",
        output:
            quant=f"{RESULTS_DIR}/nonrRNA/salmon/{{reference}}/{{sample}}/quant.sf"
        log:
            stdout=f"{RESULTS_DIR}/nonrRNA/salmon/logs/{{reference}}_{{sample}}.log"
        benchmark:
            f"{RESULTS_DIR}/nonrRNA/salmon_quant/benchmarks/{{reference}}_{{sample}}.txt"
        threads:
            config["nonrRNA"]["salmon"].get("threads", 4)
        shell:
            r"""
            set -euo pipefail

            mkdir -p $(dirname {output.quant})

            salmon quant \
                -i {input.index} \
                -l A \
                -1 {input.nonrRNA_R1} \
                -2 {input.nonrRNA_R2} \
                -p {threads} \
                -o $(dirname {output.quant}) \
                > {log.stdout} 2>&1
            """

#P4M_1  P4M_2  P4M_22  P4M_23  P4M_3  P4M_35  P4M_36  P4M_37  P4M_38  P4M_4  P4M_43  P4M_44  P4M_45  P4M_46  P4M_47  P4M_48
