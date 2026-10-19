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