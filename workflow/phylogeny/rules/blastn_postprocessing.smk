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
    wildcard_constraints:
        kingdom=BARRNAP_KINGDOM_PATTERN
    message:
        "[BLASTN] filter SILVA hits and retain top 20 references across samples reconstructed from SSU"
    input:
        blastn_tsv=f"{PHYLOGENY_DIR}/{{kingdom}}/blastn/SSU_blastn_{{kingdom}}.tsv"
    output:
        blast_top=f"{PHYLOGENY_DIR}/{{kingdom}}/blastn_postprocessing/SSU_blastn_{{kingdom}}_reference_candidates.tsv"
    params:
        top_n_hits=lambda wc: config["phylogeny"]["blastn_postprocessing"][wc.kingdom]["top_n_hits"],
        nt_identity=lambda wc: config["phylogeny"]["blastn_postprocessing"][wc.kingdom]["nt_identity"],
        query_coverage=lambda wc: config["phylogeny"]["blastn_postprocessing"][wc.kingdom]["query_coverage"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.blast_top})

        {{
            printf 'qseqid\tsseqid\tpident\tlength\tqlen\tslen\tqcovs\tevalue\tbitscore\n'
            awk '$3 >= {params.nt_identity} && $7 >= {params.query_coverage}' {input.blastn_tsv} \
                | sort -k1,1 -k9,9nr \
                | awk 'count[$1]++ < {params.top_n_hits}'
        }} > {output.blast_top}
        """


################################################################################
# CHECK WHETHER FILTERED BLAST RESULTS CONTAIN ANY HITS
################################################################################

checkpoint check_blast_top_ssu:
    wildcard_constraints:
        kingdom=BARRNAP_KINGDOM_PATTERN
    input:
        blast_top=f"{PHYLOGENY_DIR}/{{kingdom}}/blastn_postprocessing/SSU_blastn_{{kingdom}}_reference_candidates.tsv"
    output:
        status=f"{PHYLOGENY_DIR}/{{kingdom}}/blastn_postprocessing/SSU_blastn_{{kingdom}}_status.txt"
    run:
        from pathlib import Path

        with Path(input.blast_top).open(encoding="utf-8") as blast_file:
            next(blast_file, None)  # Skip the header.
            has_hits = any(line.strip() for line in blast_file)

        Path(output.status).parent.mkdir(parents=True, exist_ok=True)
        Path(output.status).write_text(
            "present\n" if has_hits else "empty\n",
            encoding="utf-8",
        )

################################################################################
# 4. Extract unique SILVA reference IDs
################################################################################

rule extract_reference_accessions:
    conda:
        "../envs/blastn.yaml"
    wildcard_constraints:
        kingdom=BARRNAP_KINGDOM_PATTERN
    message:
        "[SILVA] extract reference IDs top blast hits"
    input:
        blast_top=f"{PHYLOGENY_DIR}/{{kingdom}}/blastn_postprocessing/SSU_blastn_{{kingdom}}_reference_candidates.tsv"
    output:
        ids=f"{PHYLOGENY_DIR}/{{kingdom}}/references/samples_SSU_{{kingdom}}_accessions.txt"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.ids})

        awk -F '\t' 'NR > 1 {{print $2}}' {input.blast_top} \
            | sort -u \
            > {output.ids}
        """

################################################################################
# 5. Retrieve the reference sequences from the SILVA BLAST database
################################################################################

rule extract_reference_sequences:
    conda:
        "../envs/blastn.yaml"
    wildcard_constraints:
        kingdom=BARRNAP_KINGDOM_PATTERN
    message:
        "[SILVA] retrieve closest reference SSU sequences across samples and the reference accesions id"
    input:
        ids=f"{PHYLOGENY_DIR}/{{kingdom}}/references/samples_SSU_{{kingdom}}_accessions.txt"
    output:
        fasta=f"{PHYLOGENY_DIR}/{{kingdom}}/references/samples_SSU_{{kingdom}}_references.fasta"
    log:
        stdout=f"{PHYLOGENY_DIR}/{{kingdom}}/references/logs/blastdbcmd.log"
    benchmark:
        f"{PHYLOGENY_DIR}/{{kingdom}}/references/benchmarks/blastdbcmd.txt"
    params:
        db=lambda wc: config["phylogeny"]["blastn"]["databases"][wc.kingdom],
    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {output.fasta})

        blastdbcmd \
            -db {params.db} \
            -entry_batch {input.ids} \
            -outfmt "%f" \
            > {output.fasta} \
            2> {log.stdout}
        """

################################################################################
# 6. Rename SILVA reference sequences
################################################################################

rule rename_reference_headers:
    wildcard_constraints:
        kingdom=BARRNAP_KINGDOM_PATTERN
    message:
        "Keep only reference accessions in FASTA headers"
    input:
        orig_fasta=f"{PHYLOGENY_DIR}/{{kingdom}}/references/samples_SSU_{{kingdom}}_references.fasta"
    output:
        renamed_fasta=f"{PHYLOGENY_DIR}/{{kingdom}}/references/samples_SSU_{{kingdom}}_renamed.fasta"
    params:
        lineages=lambda wc: config["phylogeny"]["blastn"]["lineage"][wc.kingdom],
    shell:
        r"""
        set -euo pipefail

        python {PHYLOGENEY_SCRIPTS_DIR}/reorder_header.py --in_fasta {input.orig_fasta} --taxonomy_tsv {params.lineages} --out_fasta {output.renamed_fasta}
        """

################################################################################
# 6. COMBINE RECONSTRUCTED SSUs AND SILVA REFERENCES
################################################################################

rule combine_ssu_references:
    wildcard_constraints:
        kingdom=BARRNAP_KINGDOM_PATTERN
    message:
        "[Phylogeny] Combine reconstructed SSUs and SILVA references"
    input:
        query=f"{ASSEMBLY_DIR}/SSU_candidates/postprocessing/{{kingdom}}/SSU_{{kingdom}}_catalogue.fasta",
        references=f"{PHYLOGENY_DIR}/{{kingdom}}/references/samples_SSU_{{kingdom}}_references.fasta"
    output:
        fasta=f"{PHYLOGENY_DIR}/{{kingdom}}/references/samples_SSU_{{kingdom}}_with_references.fasta"
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname "{output.fasta}")"

        cat \
            "{input.query}" \
            "{input.references}" \
            > "{output.fasta}"
        """