################################################################################
# FILTER BLAST RESULTS AND RETAIN TOP REFERENCES PER QUERY
################################################################################

rule filter_blast_ssu:
    conda:
        "../envs/blastn.yaml"
    wildcard_constraints:
        domain=PHYLOGENY_DOMAIN_PATTERN
    message:
        "[BLASTN] Selecting reference candidates for {wildcards.domain}"
    input:
        blastn_tsv=f"{PHYLOGENY_DIR}/{{domain}}/blastn/SSU_blastn_{{domain}}.tsv"
    output:
        blast_top=f"{PHYLOGENY_DIR}/{{domain}}/blastn_postprocessing/SSU_blastn_{{domain}}_reference_candidates.tsv"
    params:
        top_n_hits=lambda wc: config["phylogeny"]["blastn_postprocessing"][wc.domain]["top_n_hits"],
        nt_identity=lambda wc: config["phylogeny"]["blastn_postprocessing"][wc.domain]["nt_identity"],
        query_coverage=lambda wc: config["phylogeny"]["blastn_postprocessing"][wc.domain]["query_coverage"]
    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {output.blast_top})

        {{
            printf 'qseqid\tsseqid\tpident\tlength\tqlen\tslen\tqcovs\tevalue\tbitscore\n'

            awk '$3 >= {params.nt_identity} && $7 >= {params.query_coverage}' {input.blastn_tsv} \
                | LC_ALL=C sort -k1,1 -k9,9nr -k2,2 \
                | awk 'count[$1]++ < {params.top_n_hits}'
        }} > {output.blast_top}
        """


################################################################################
# CHECK WHETHER FILTERED BLAST RESULTS CONTAIN HITS
################################################################################

checkpoint check_blast_top_ssu:
    wildcard_constraints:
        domain=PHYLOGENY_DOMAIN_PATTERN
    input:
        blast_top=f"{PHYLOGENY_DIR}/{{domain}}/blastn_postprocessing/SSU_blastn_{{domain}}_reference_candidates.tsv"
    output:
        status=f"{PHYLOGENY_DIR}/{{domain}}/blastn_postprocessing/SSU_blastn_{{domain}}_status.txt"
    run:
        with open(input.blast_top) as handle:
            next(handle, None)
            has_hits = any(line.strip() for line in handle)

        Path(output.status).parent.mkdir(parents=True, exist_ok=True)
        Path(output.status).write_text(
            "present\n" if has_hits else "empty\n",
            encoding="utf-8",
        )


################################################################################
# EXTRACT UNIQUE REFERENCE IDS
################################################################################

rule extract_reference_accessions:
    conda:
        "../envs/blastn.yaml"
    wildcard_constraints:
        domain=PHYLOGENY_DOMAIN_PATTERN
    message:
        "[References] Extracting reference IDs for {wildcards.domain}"
    input:
        blast_top=f"{PHYLOGENY_DIR}/{{domain}}/blastn_postprocessing/SSU_blastn_{{domain}}_reference_candidates.tsv"
    output:
        ids=f"{PHYLOGENY_DIR}/{{domain}}/references/samples_SSU_{{domain}}_accessions.txt"
    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {output.ids})

        awk -F '\t' 'NR > 1 {{print $2}}' {input.blast_top} \
            | LC_ALL=C sort -u \
            > {output.ids}
        """


################################################################################
# RETRIEVE REFERENCE SEQUENCES
################################################################################

rule extract_reference_sequences:
    conda:
        "../envs/blastn.yaml"
    wildcard_constraints:
        domain=PHYLOGENY_DOMAIN_PATTERN
    message:
        "[References] Retrieving SSU references for {wildcards.domain}"
    input:
        ids=f"{PHYLOGENY_DIR}/{{domain}}/references/samples_SSU_{{domain}}_accessions.txt"
    output:
        fasta=f"{PHYLOGENY_DIR}/{{domain}}/references/samples_SSU_{{domain}}_references.fasta"
    log:
        stdout=f"{PHYLOGENY_DIR}/{{domain}}/references/logs/blastdbcmd.log"
    benchmark:
        f"{PHYLOGENY_DIR}/{{domain}}/references/benchmarks/blastdbcmd.txt"
    params:
        db=lambda wc: config["phylogeny"]["blastn"]["databases"][wc.domain]
    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {output.fasta}) $(dirname {log.stdout})

        blastdbcmd \
            -db {params.db} \
            -entry_batch {input.ids} \
            -outfmt "%f" \
            > {output.fasta} \
            2> {log.stdout}
        """


################################################################################
# RENAME REFERENCE HEADERS
################################################################################

rule rename_reference_headers:
    wildcard_constraints:
        domain=PHYLOGENY_DOMAIN_PATTERN
    message:
        "[References] Renaming reference headers for {wildcards.domain}"
    input:
        orig_fasta=f"{PHYLOGENY_DIR}/{{domain}}/references/samples_SSU_{{domain}}_references.fasta"
    output:
        renamed_fasta=f"{PHYLOGENY_DIR}/{{domain}}/references/samples_SSU_{{domain}}_renamed.fasta"
    params:
        script=f"{PHYLOGENY_SCRIPTS_DIR}/reorder_header.py",
        lineages=lambda wc: config["phylogeny"]["blastn"]["lineage"][wc.domain]
    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {output.renamed_fasta})

        python3 {params.script} \
            --in_fasta {input.orig_fasta} \
            --taxonomy_tsv {params.lineages} \
            --out_fasta {output.renamed_fasta}
        """


################################################################################
# COMBINE SAMPLE-DERIVED SSUS AND RENAMED REFERENCES
################################################################################

rule combine_ssu_references:
    wildcard_constraints:
        domain=PHYLOGENY_DOMAIN_PATTERN
    message:
        "[Phylogeny] Combining SSUs and references for {wildcards.domain}"
    input:
        query=f"{ASSEMBLY_DIR}/cross_assembly/SSU/domains/{{domain}}/final/transcript.fasta",
        references=f"{PHYLOGENY_DIR}/{{domain}}/references/samples_SSU_{{domain}}_renamed.fasta"
    output:
        fasta=f"{PHYLOGENY_DIR}/{{domain}}/references/samples_SSU_{{domain}}_with_references.fasta"
    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {output.fasta})

        cat \
            {input.query} \
            {input.references} \
            > {output.fasta}
        """