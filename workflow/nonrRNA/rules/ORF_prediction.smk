rule Prodigal:
    conda:
        "../envs/prodigal.yaml"
    message:
        "[Prodigal] ORF prediction of microbial (bacterial and archaeal) genes"
    input:
        assembled_fq=f"{RESULTS_DIR}/nonrRNA/assembly/nonrRNA_assembly.fasta"
    output:
        predicted_genes=f"{RESULTS_DIR}/nonrRNA/predicted/ORF_genes.gff",
        predicted_genes_nt=f"{RESULTS_DIR}/nonrRNA/predicted/ORF_genes.fasta"
        predicted_proteins=f"{RESULTS_DIR}/nonrRNA/predicted/ORF_proteins.faa"
    log:
        f"{RESULTS_DIR}/nonrRNA/predicted/logs/prodigal_predicted.log"
    benchmark:
        f"{RESULTS_DIR}/nonrRNA/predicted/benchmarks/prodigal_predicted.txt"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.predicted_genes})

        prodigal -i {input.assembled_fq} \
            -a {output.predicted_proteins} \
            -d {output.predicted_genes_nt} \
            -o {output.predicted_ORF} \
            -p meta -f gff -m > {log} 2>&1
        """

# Transdecoder are usually best for eukaryotes

#extract candidates
rule TransDecoderLongOrfs:
    conda:
        "../envs/transdecoder.yaml"
    message:
        "[TransDecoder.LongOrfs] Identify candidate open reading frames"
    input:
        assembled_fq=f"{RESULTS_DIR}/nonrRNA/assembly/nonrRNA_assembly.fasta"
    output:
        gene_map=f"{RESULTS_DIR}/nonrRNA/predicted/gene_to_transcripts.tsv",
    log:
        f"{RESULTS_DIR}/nonrRNA/predicted/logs/prodigal_predicted.log"
    benchmark:
        f"{RESULTS_DIR}/nonrRNA/predicted/benchmarks/prodigal_predicted.txt"
    params:
        result_dir=f"{RESULTS_DIR}/nonrRNA/predicted/",
        options="-m 100"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.gene_map})

        TransDecoder.LongOrfs -t {input.assembled_fq} \
            --gene_trans_map {output.gene_map} \
            --output_dir {params.result_dir} \
            {params.options}
        """

#Score ORFs to select most likely true coding sequences 
rule TransDecoderPredict:
    conda:
        "../envs/transdecoder.yaml"
    message:
        "[TransDecoder.Predict] score ORF candidates and select most likely true coding sequence"
    input:
        assembled_fq=f"{RESULTS_DIR}/nonrRNA/assembly/nonrRNA_assembly.fasta",
    output:
        pfam_hits=f"{RESULTS_DIR}/nonrRNA/predicted/transdecoder_pfam.tsv"
        blastp_hits=f"{RESULTS_DIR}/nonrRNA/predicted/transdecoder_blastp.tsv"
    log:
        f"{RESULTS_DIR}/nonrRNA/predicted/logs/transdecoder_predicted.log"
    benchmark:
        f"{RESULTS_DIR}/nonrRNA/predicted/benchmarks/transdecoder_predicted.txt"
    params:
        result_dir=f"{RESULTS_DIR}/nonrRNA/predicted/",
        options="-m 100"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.pfam_hits})

        TransDecoder.LongOrfs -t {input.assembled_fq} \
            --output_dir {params.result_dir} \
            --retain_pfam_hits {output.pfam_hits} \
            --retain_blastp_hits {output.blastp_hits}
        """

# TransDecoder.LongOrfs -t "contigs.fasta"
# TransDecoder.Predict -t "contigs.fasta"