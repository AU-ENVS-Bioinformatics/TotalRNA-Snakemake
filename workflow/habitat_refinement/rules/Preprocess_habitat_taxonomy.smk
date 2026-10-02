rule preprocess_habitat_taxonomy:
    conda:
        str(TAXONOMY_ENVS_DIR / "krakentools.yaml")
    message:
        "[KrakenTools] Extracting shared glacier-order reads from {wildcards.sample}"
    input:
        read_1=lambda wc: grouped_ssu_read(wc.sample, "prokaryotes", 1),
        read_2=lambda wc: grouped_ssu_read(wc.sample, "prokaryotes", 2),
        report=f"{TAXONOMY_DIR}/{{sample}}/SSU/silva/kraken/{{sample}}.silva.k2report",
        kraken=f"{TAXONOMY_DIR}/{{sample}}/SSU/silva/kraken/{{sample}}.silva.kraken",
        taxids=config["Habitat_refinement"]["Glacier"]["order_taxids"],
    output:
        read_1=f"{HABITAT_ANALYSIS_DIR}/{{sample}}/Glacier/grouped_reads/{{sample}}.glacial_orders_R1.fastq.gz",
        read_2=f"{HABITAT_ANALYSIS_DIR}/{{sample}}/Glacier/grouped_reads/{{sample}}.glacial_orders_R2.fastq.gz",
    log:
        stdout=f"{HABITAT_ANALYSIS_DIR}/{{sample}}/Glacier/logs/{{sample}}.preprocess_habitat_taxonomy.log",
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname {output.read_1})"

        (
            # Reject empty files and nonnumeric taxids.
            awk '
                {{
                    sub(/\r$/, "")
                    if ($0 !~ /^[0-9]+$/) {{
                        printf "Invalid taxid at line %d: %s\n", NR, $0 \
                            > "/dev/stderr"
                        invalid = 1
                    }}
                }}
                END {{
                    if (NR == 0) {{
                        print "Taxid file is empty" > "/dev/stderr"
                        invalid = 1
                    }}
                    exit invalid
                }}
            ' {input.taxids:q}

            mapfile -t taxids < <(
                sed 's/\r$//' {input.taxids:q} | LC_ALL=C sort -u
            )

            extract_kraken_reads.py \
                -k {input.kraken} \
                -r {input.report} \
                -s1 {input.read_1} \
                -s2 {input.read_2} \
                -t "${{taxids[@]}}" \
                --include-children \
                --fastq-output \
                --max 10000000000 \
                -o {output.read_1}.tmp.fastq \
                -o2 {output.read_2}.tmp.fastq

            gzip -c {output.read_1}.tmp.fastq > {output.read_1}
            gzip -c {output.read_2}.tmp.fastq > {output.read_2}

            rm -- \
                {output.read_1}.tmp.fastq \
                {output.read_2}.tmp.fastq
        ) > {log.stdout} 2>&1
        """