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