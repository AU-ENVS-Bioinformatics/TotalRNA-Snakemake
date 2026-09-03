################################################################################
# ITS BLASTN TAXONOMY AGAINST UNITE
################################################################################

if its_taxonomy_enabled(config):

    ITS_REGIONS = its_blast_regions(config)
    ITS_REGION_PATTERN = "|".join(ITS_REGIONS)

    rule blastn_its_unite:
        conda:
            "../envs/blastn.yaml"
        wildcard_constraints:
            region=ITS_REGION_PATTERN
        message:
            "[BLASTN] Classifying {wildcards.region} sequences from {wildcards.sample} against UNITE"
        input:
            query=f"{ASSEMBLY_DIR}/{{sample}}/ITS_candidates/ITSx/{{sample}}_ITSx.{{region}}.fasta"
        output:
            hits=f"{TAXONOMY_DIR}/{{sample}}/ITS/{{region}}/blastn/{{sample}}.{{region}}.UNITE.blastn.tsv"
        log:
            stdout=f"{TAXONOMY_DIR}/{{sample}}/ITS/{{region}}/blastn/logs/{{sample}}.{{region}}.UNITE.blastn.log"
        benchmark:
            f"{TAXONOMY_DIR}/{{sample}}/ITS/{{region}}/blastn/benchmarks/{{sample}}.{{region}}.UNITE.blastn.txt"
        params:
            db=config["taxonomy"]["ITS"]["blastn"]["database"],
            options=config["taxonomy"]["ITS"]["blastn"].get("options", ""),
        threads:
            config["taxonomy"]["ITS"]["blastn"].get("threads", 12)
        shell:
            r"""
            set -euo pipefail

            mkdir -p $(dirname {output.hits})

            if [[ -s {input.query} ]]; then
                blastn \
                    -query {input.query} \
                    -db {params.db} \
                    -num_threads {threads} \
                    -out {output.hits} \
                    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen qcovs qcovhsp evalue bitscore stitle" \
                    {params.options} \
                    > {log.stdout} 2>&1
            else
                touch {output.hits}

                echo \
                    "[BLASTN] Empty input FASTA: {input.query}" \
                    > {log.stdout}
            fi
            """