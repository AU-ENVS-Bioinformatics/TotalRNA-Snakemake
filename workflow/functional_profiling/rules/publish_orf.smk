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