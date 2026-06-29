PREDICTED_DIR = f"{RESULTS_DIR}/nonrRNA/predicted"

if config["nonrRNA"]["module"] == "coassembly" and config["nonrRNA"]["predictor"] == "prodigal":
    rule Prodigal:
        conda:
            "../envs/prodigal.yaml"
        message:
            "[Prodigal] ORF prediction of microbial (bacterial and archaeal) genes"
        input:
            assembled_fa=f"{RESULTS_DIR}/nonrRNA/coassembly/transcripts.fasta"
        output:
            predicted_genes_nt=f"{PREDICTED_DIR}/ORF_genes.fasta",
            predicted_proteins=f"{PREDICTED_DIR}/ORF_proteins.faa",
            predicted_genes=f"{PREDICTED_DIR}/ORF_genes.gff",
            predicted_genes_scores=f"{PREDICTED_DIR}/ORF_genes_scores.tsv"
        log:
            stdout = f"{PREDICTED_DIR}/logs/prodigal_predicted.log"
        benchmark:
            f"{PREDICTED_DIR}/benchmarks/prodigal_predicted.txt"
        shell:
            r"""
            set -euo pipefail
            mkdir -p $(dirname {output.predicted_genes})

            prodigal -i {input.assembled_fa} \
                -a {output.predicted_proteins} \
                -d {output.predicted_genes_nt} \
                -o {output.predicted_ORF} \
                -s {output.predicted_genes_scoresscre} \
                > {log.stdout} 2>&1
            """
elif config["nonrRNA"]["module"] == "coassembly" and config["nonrRNA"]["predictor"] == "transdecoder":
    # Transdecoder are usually best for eukaryotes
    
    rule TransDecoderLongOrfs:
        conda:
            "../envs/transdecoder.yaml"
        message:
            "[TransDecoder.LongOrfs] Identify candidate ORFs"
        input:
            assembled_fa=f"{RESULTS_DIR}/nonrRNA/coassembly/transcripts.fasta"
        output:
            pep=temp(f"{PREDICTED_DIR}/longest_orfs.pep"),
            pep2=f"{PREDICTED_DIR}/longest_orfs_tmp.pep",
            cds=temp(f"{PREDICTED_DIR}/longest_orfs.cds"),
            gff=temp(f"{PREDICTED_DIR}/longest_orfs.gff3"),
        log:
            stdout = f"{PREDICTED_DIR}/logs/longorfs.log"
        benchmark:
            f"{RESULTS_DIR}/nonrRNA/predicted/benchmarks/longorfs.txt"
        params:
            outdir=f"{PREDICTED_DIR}",
            options=config["nonrRNA"]["transdecoder"]["options"],
        shell:
            r"""
            set -euo pipefail
            mkdir -p $(dirname {output.pep})

            TransDecoder.LongOrfs \
                -t {input.assembled_fa} \
                --output_dir {params.outdir} \
                {params.options} \
                > {log.stdout} 2>&1
            

            TD_DIR="{params.outdir}/$(basename {input.assembled_fa}).transdecoder_dir"
            
            ln -sf "$TD_DIR/longest_orfs.pep" {output.pep}
            ln -sf "$TD_DIR/longest_orfs.cds" {output.cds}
            ln -sf "$TD_DIR/longest_orfs.gff3" {output.gff}

            cp "$TD_DIR/longest_orfs.pep" {output.pep2}
            """

    # TransDecoder.LongOrfs -t transcripts.fasta --output_dir ../Transdecoder
    #-m 100 \--single_best_only 
    
    #Score ORFs to select most likely true coding sequences 
    rule TransDecoderPredict:
        conda:
            "../envs/transdecoder.yaml"
        message:
            "[TransDecoder.Predict] score ORF candidates and select most likely true coding sequence"
        input:
            assembled_fa=f"{RESULTS_DIR}/nonrRNA/coassembly/transcripts.fasta",
            longorfs=f"{RESULTS_DIR}/nonrRNA/predicted/longest_orfs.pep",
            diamond=f"{RESULTS_DIR}/nonrRNA/diamond/transdecoder.blastp.outfmt6",
            pfam=f"{RESULTS_DIR}/nonrRNA/hmmscan/pfam.domtblout",
        output:
            pep = f"{PREDICTED_DIR}/transcripts.fasta.transdecoder.pep",
            cds = f"{PREDICTED_DIR}/transcripts.fasta.transdecoder.cds",
            gff = f"{PREDICTED_DIR}/transcripts.fasta.transdecoder.gff3",
            done=f"{PREDICTED_DIR}/transdecoder.done"
        log:
            stdout = f"{RESULTS_DIR}/nonrRNA/predicted/logs/transdecoder_predicted.log"
        benchmark:
            f"{RESULTS_DIR}/nonrRNA/predicted/benchmarks/transdecoder_predicted.txt"
        params:
            outdir=f"{PREDICTED_DIR}",
        shell:
            r"""
            set -euo pipefail
            mkdir -p $(dirname {output.done})

            TransDecoder.Predict \
                -t {input.assembled_fa} \
                --retain_blastp_hits {input.diamond} \
                --retain_pfam_hits {input.pfam} \
                --single_best_only \
                --output_dir {params.outdir} \
                > {log} 2>&1

            touch {output.done}
            """

    # TransDecoder.LongOrfs -t "contigs.fasta"
    # TransDecoder.Predict -t "contigs.fasta"

# filtering ORFs: keep biologically meaningful coding sequences 
# using:
# homology evidence (DIAMOND)
# domain presence (Pfam)