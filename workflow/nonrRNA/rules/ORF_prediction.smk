################################################################################
# ORF PREDICTION
################################################################################

FUNCTIONAL_MODULE = config["functional_profiling"]["module"]
ORF_PREDICTOR = config["functional_profiling"]["ORF_predictor"]

ASSEMBLED_NONRRNA = f"{ASSEMBLY_DIR}/coassembly/non_rRNA/rnaspades/transcripts.fasta"

ORF_DIR = f"{FUNCTION_DIR}/ORF_prediction"
TRANSDECODER_WORK_DIR = f"{ORF_DIR}/transdecoder/work"
TRANSDECODER_LONGORFS_DIR = f"{ORF_DIR}/transdecoder/longorfs"
TRANSDECODER_PREDICT_DIR = f"{ORF_DIR}/transdecoder/predict"

SUPPORTED_ORF_PREDICTORS = {
    "prodigal",
    "transdecoder",
}

if ORF_PREDICTOR not in SUPPORTED_ORF_PREDICTORS:
    raise ValueError(
        f"Unsupported functional_profiling.ORF_predictor: {ORF_PREDICTOR}. "
        f"Supported values are: {', '.join(sorted(SUPPORTED_ORF_PREDICTORS))}"
    )


################################################################################
# PRODIGAL
################################################################################

if FUNCTIONAL_MODULE in ["assembly", "both", "all"]:

    if ORF_PREDICTOR == "prodigal":

        rule prodigal_orf_prediction:
            conda:
                "../envs/prodigal.yaml"
            message:
                "[Prodigal] Predict bacterial and archaeal ORFs from the non-rRNA coassembly"
            input:
                transcripts=ASSEMBLED_NONRRNA
            output:
                proteins=f"{ORF_DIR}/prodigal/ORF_proteins.faa",
                genes=f"{ORF_DIR}/prodigal/ORF_genes.fasta",
                gff=f"{ORF_DIR}/prodigal/ORF_genes.gff3",
                scores=f"{ORF_DIR}/prodigal/ORF_genes_scores.sco"
            log:
                stdout=f"{ORF_DIR}/prodigal/logs/prodigal.log"
            benchmark:
                f"{ORF_DIR}/prodigal/benchmarks/prodigal.txt"
            params:
                options=config["functional_profiling"]["prodigal"].get("options", "")
            shell:
                r"""
                set -euo pipefail

                mkdir -p $(dirname {output.proteins})

                prodigal \
                    -i {input.transcripts} \
                    -a {output.proteins} \
                    -d {output.genes} \
                    -o {output.gff} \
                    -s {output.scores} \
                    {params.options} \
                    > {log.stdout} 2>&1
                """


################################################################################
# TRANSDECODER.LONGORFS
################################################################################

    elif ORF_PREDICTOR == "transdecoder":

        rule transdecoder_longorfs:
            conda:
                "../envs/transdecoder.yaml"

            message:
                "[TransDecoder.LongOrfs] Identify candidate ORFs in the non-rRNA coassembly"

            input:
                transcripts=ASSEMBLED_NONRRNA

            output:
                pep=f"{TRANSDECODER_LONGORFS_DIR}/longest_orfs.pep",
                cds=f"{TRANSDECODER_LONGORFS_DIR}/longest_orfs.cds",
                gff=f"{TRANSDECODER_LONGORFS_DIR}/longest_orfs.gff3"

            log:
                stdout=f"{TRANSDECODER_LONGORFS_DIR}/logs/longorfs.log"

            benchmark:
                f"{TRANSDECODER_LONGORFS_DIR}/benchmarks/longorfs.txt"

            params:
                workdir=TRANSDECODER_WORK_DIR,
                options=config["functional_profiling"]["transdecoder"].get("options", "")

            shell:
                r"""
                set -euo pipefail

                mkdir -p \
                    {params.workdir} \
                    $(dirname {output.pep})

                TransDecoder.LongOrfs \
                    -t {input.transcripts} \
                    --output_dir {params.workdir} \
                    {params.options} \
                    > {log.stdout} 2>&1

                TD_DIR={params.workdir}/$(basename {input.transcripts}).transdecoder_dir

                ln -sfn $(realpath $TD_DIR/longest_orfs.pep) {output.pep}
                ln -sfn $(realpath $TD_DIR/longest_orfs.cds) {output.cds}
                ln -sfn $(realpath $TD_DIR/longest_orfs.gff3) {output.gff}
                """


################################################################################
# TRANSDECODER.PREDICT
################################################################################

        rule transdecoder_predict:
            conda:
                "../envs/transdecoder.yaml"
            message:
                "[TransDecoder.Predict] Select ORFs using DIAMOND and Pfam evidence"
            input:
                transcripts=ASSEMBLED_NONRRNA,
                longorfs=f"{TRANSDECODER_LONGORFS_DIR}/longest_orfs.pep",
                diamond=f"{FUNCTION_DIR}/diamond/transdecoder/transdecoder.blastp.outfmt6",
                pfam=f"{FUNCTION_DIR}/hmmscan/merged/pfam.domtblout"

            output:
                pep=f"{TRANSDECODER_PREDICT_DIR}/transcripts.fasta.transdecoder.pep",
                cds=f"{TRANSDECODER_PREDICT_DIR}/transcripts.fasta.transdecoder.cds",
                gff=f"{TRANSDECODER_PREDICT_DIR}/transcripts.fasta.transdecoder.gff3"

            log:
                stdout=f"{TRANSDECODER_PREDICT_DIR}/logs/predict.log"

            benchmark:
                f"{TRANSDECODER_PREDICT_DIR}/benchmarks/predict.txt"

            params:
                workdir=TRANSDECODER_WORK_DIR

            shell:
                r"""
                set -euo pipefail

                mkdir -p \
                    {params.workdir} \
                    $(dirname {output.pep})

                TransDecoder.Predict \
                    -t {input.transcripts} \
                    --retain_blastp_hits {input.diamond} \
                    --retain_pfam_hits {input.pfam} \
                    --single_best_only \
                    --output_dir {params.workdir} \
                    > {log.stdout} 2>&1

                ln -sfn $(realpath {params.workdir}/transcripts.fasta.transdecoder.pep) {output.pep}
                ln -sfn $(realpath {params.workdir}/transcripts.fasta.transdecoder.cds) {output.cds}
                ln -sfn $(realpath {params.workdir}/transcripts.fasta.transdecoder.gff3) {output.gff}
                """


################################################################################
# STANDARDIZED FINAL ORF OUTPUTS
################################################################################

if FUNCTIONAL_MODULE in ["assembly", "both", "all"]:

    if ORF_PREDICTOR == "prodigal":
        FINAL_ORF_PROTEINS = f"{ORF_DIR}/prodigal/ORF_proteins.faa"
        FINAL_ORF_GENES = f"{ORF_DIR}/prodigal/ORF_genes.fasta"
        FINAL_ORF_GFF = f"{ORF_DIR}/prodigal/ORF_genes.gff3"

    elif ORF_PREDICTOR == "transdecoder":
        FINAL_ORF_PROTEINS = f"{TRANSDECODER_PREDICT_DIR}/transcripts.fasta.transdecoder.pep"
        FINAL_ORF_GENES = f"{TRANSDECODER_PREDICT_DIR}/transcripts.fasta.transdecoder.cds"
        FINAL_ORF_GFF = f"{TRANSDECODER_PREDICT_DIR}/transcripts.fasta.transdecoder.gff3"


    rule publish_final_orfs:
        message:
            f"[ORF prediction] Publish standardized {ORF_PREDICTOR} outputs"
        input:
            proteins=FINAL_ORF_PROTEINS,
            genes=FINAL_ORF_GENES,
            gff=FINAL_ORF_GFF
        output:
            proteins=f"{ORF_DIR}/final_orf/ORF_proteins.faa",
            genes=f"{ORF_DIR}/final_orf/ORF_genes.fasta",
            gff=f"{ORF_DIR}/final_orf/ORF_genes.gff3"
        shell:
            r"""
            set -euo pipefail

            mkdir -p $(dirname {output.proteins})

            ln -sfn $(realpath {input.proteins}) {output.proteins}
            ln -sfn $(realpath {input.genes}) {output.genes}
            ln -sfn $(realpath {input.gff}) {output.gff}
            """