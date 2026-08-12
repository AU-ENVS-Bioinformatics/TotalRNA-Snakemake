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

rule filter_blast_ssu:
    conda:
        "../envs/blastn.yaml"
    message:
        "[BLASTN] filter SILVA hits and retain top 20 references per SSU for {wildcards.sample}"
    input:
        blast=f"{PHYLOGENY_DIR}/{{sample}}/blastn/{{sample}}_SSU_blastn.tsv"
    output:
        blast_top=f"{PHYLOGENY_DIR}/{{sample}}/blastn/{{sample}}_SSU_blast_top20.tsv"
    params:
        top_n_hits=config["phylogeny"]["blastn_postprocessing"]["top_n_hits"],
        nt_identify=config["phylogeny"]["blastn_postprocessing"]["nt_identify"],
        query_coverage=config["phylogeny"]["blastn_postprocessing"]["query_coverage"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.blast_top})

        awk '$3 >= {params.nt_identify} && $7 >= {params.query_coverage}' {input.blast} | sort -k1,1 -k9,9nr | awk 'count[$1]++ < {params.top_n_hits}' > {output.blast_top}
        """

################################################################################
# 4. Extract unique SILVA reference IDs
################################################################################

rule extract_silva_accessions:
    conda:
        "../envs/blastn.yaml"
    message:
        "[SILVA] extract reference IDs for {wildcards.sample}"
    input:
        filtered=f"{PHYLOGENY_DIR}/{{sample}}/blastn/{{sample}}_SSU_blast_top20.tsv"
    output:
        ids=f"{PHYLOGENY_DIR}/{{sample}}/references/{{sample}}_SILVA_accessions.txt"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.ids})

        cut -f2 {input.filtered} \
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
        "[SILVA] retrieve closest reference SSU sequences for {wildcards.sample}"
    input:
        ids=f"{PHYLOGENY_DIR}/{{sample}}/references/{{sample}}_SILVA_accessions.txt"
    output:
        fasta=f"{PHYLOGENY_DIR}/{{sample}}/references/{{sample}}_SILVA_references.fasta"
    log:
        f"{PHYLOGENY_DIR}/{{sample}}/logs/{{sample}}_blastdbcmd.log"
    benchmark:
        f"{PHYLOGENY_DIR}/{{sample}}/benchmarks/{{sample}}_blastdbcmd.txt"
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
# 6. Combine phyloFlash sequences and SILVA reference sequences
################################################################################

rule combine_ssu_references:
    conda:
        "../envs/blastn.yaml"
    message:
        "[Phylogeny] combine reconstructed SSUs and SILVA references for {wildcards.sample}"
    input:
        query=f"{PHYLOGENY_DIR}/{{sample}}/phyloflash/{{sample}}.all.final.fasta",
        refs=f"{PHYLOGENY_DIR}/{{sample}}/references/{{sample}}_SILVA_references.fasta"
    output:
        fasta=f"{PHYLOGENY_DIR}/{{sample}}/phyloflashreference/{{sample}}_SSU_with_references.fasta"
    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {output.fasta})

        cat \
            {input.query} \
            {input.refs} \
            > {output.fasta}
        """
