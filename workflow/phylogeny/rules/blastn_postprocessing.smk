################################################################################
# 3. Filter BLAST results
#
# Requirements:
#     >= 97% nucleotide identity
#     >= 90% query coverage
#
# Then retain the best 20 hits PER reconstructed SSU sequence.
#
# Columns:
# 1 qseqid
# 2 sseqid
# 3 pident
# 4 length
# 5 qlen
# 6 slen
# 7 qcovs
# 8 evalue
# 9 bitscore
################################################################################

PHYLOGENEY_SCRIPTS_DIR = str(Path(workflow.basedir)/"phylogeny/scripts")

rule filter_blast_ssu:
    conda:
        "../envs/blastn.yaml"
    message:
        "[BLASTN] filter SILVA hits and retain top 20 references across samples reconstructed from SSU"
    input:
        blastn_tsv=f"{PHYLOGENY_DIR}/cross_sample/blastn/cross_sample_SSU_blastn.tsv"
    output:
        blast_top=f"{PHYLOGENY_DIR}/cross_sample/blastn/cross_sample_SSU_blast_top20.tsv"
    params:
        top_n_hits=config["phylogeny"]["blastn_postprocessing"]["top_n_hits"],
        nt_identify=config["phylogeny"]["blastn_postprocessing"]["nt_identify"],
        query_coverage=config["phylogeny"]["blastn_postprocessing"]["query_coverage"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.blast_top})

        awk '$3 >= {params.nt_identify} && $7 >= {params.query_coverage}' {input.blastn_tsv} | sort -k1,1 -k9,9nr | awk 'count[$1]++ < {params.top_n_hits}' > {output.blast_top}
        """

################################################################################
# 4. Extract unique SILVA reference IDs
################################################################################

rule extract_silva_accessions:
    conda:
        "../envs/blastn.yaml"
    message:
        "[SILVA] extract reference IDs top blast hits"
    input:
        blast_top=f"{PHYLOGENY_DIR}/cross_sample/blastn/cross_sample_SSU_blast_top20.tsv"
    output:
        ids=f"{PHYLOGENY_DIR}/references/cross_sample_SSU_SILVA_accessions.txt"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.ids})

        cut -f2 {input.blast_top} \
            | sort -u \
            > {output.ids}
        """

################################################################################
# 5. Retrieve the reference sequences from the SILVA BLAST database
################################################################################

rule extract_silva_sequences:
    conda:
        "../envs/blastn.yaml"
    message:
        "[SILVA] retrieve closest reference SSU sequences across samples and the reference accesions id"
    input:
        ids=f"{PHYLOGENY_DIR}/references/cross_sample_SSU_SILVA_accessions.txt"
    output:
        fasta=f"{PHYLOGENY_DIR}/references/cross_sample_acc_SSU_SILVA_references.fasta"
    log:
        f"{PHYLOGENY_DIR}/cross_sample/logs/blastdbcmd.log"
    benchmark:
        f"{PHYLOGENY_DIR}/cross_sample/benchmarks/blastdbcmd.txt"
    params:
        db=config["databases"]["SILVA_138.2_N99_blast"],
    shell:
        r"""
        set -euo pipefail

        blastdbcmd \
            -db {params.db} \
            -entry_batch {input.ids} \
            -outfmt "%f" \
            > {output.fasta} \
            2> {log}
        """

################################################################################
# 6. Rename SILVA reference sequences
################################################################################

rule rename_silva_references:
    message:
        "Reorder SILVA reference headers to accomodate taxonomical ranking"
    input:
        orig_fasta=f"{PHYLOGENY_DIR}/references/cross_sample_acc_SSU_SILVA_references.fasta"
    output:
        renamed_fasta=f"{PHYLOGENY_DIR}/references/cross_sample_acc_SSU_renamed.fasta"
    shell:
        r"""
        set -euo pipefail

        python {PHYLOGENEY_SCRIPTS_DIR}/reorder_header.py --in_fasta {input.orig_fasta} --out_fasta {output.renamed_fasta}
        """

################################################################################
# 6. Combine phyloFlash sequences and SILVA reference sequences
################################################################################

rule combine_ssu_references:
    conda:
        "../envs/blastn.yaml"
    message:
        "[Phylogeny] combine reconstructed SSUs and SILVA references across samples"
    input:
        fasta_query=f"{PHYLOGENY_DIR}/cross_sample/reconstructed_SSU_all_samples.fasta",
        fasta_ref=f"{PHYLOGENY_DIR}/references/cross_sample_acc_SSU_renamed.fasta"
    output:
        fasta=f"{PHYLOGENY_DIR}/phyloflashreference/cross_sample_acc_SSU_with_references.fasta"
    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {output.fasta})

        cat \
            {input.fasta_query} \
            {input.fasta_ref} \
            > {output.fasta}
        """
