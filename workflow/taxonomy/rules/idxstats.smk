rule taxonomy_bwamem2_idxstats:
    conda:
        "../envs/bwamem2_samtools.yaml"
    wildcard_constraints:
        marker=TAXONOMY_ASSEMBLY_MARKER_PATTERN
    input:
        bam=rules.taxonomy_bwamem2_mem.output.bam
    output:
        stats=f"{TAXONOMY_ASSEMBLY_DIR}/{{sample}}/{{marker}}/alignment/{{sample}}.idxstats.tsv"
    log:
        stderr=f"{TAXONOMY_ASSEMBLY_DIR}/{{sample}}/{{marker}}/logs/{{sample}}.idxstats.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.stats}")" "$(dirname "{log.stderr}")"
        samtools idxstats "{input.bam}" > "{output.stats}" 2> "{log.stderr}"
        """

rule taxonomy_contig_count_matrix:
    wildcard_constraints:
        marker=TAXONOMY_ASSEMBLY_MARKER_PATTERN
    input:
        stats=lambda wc: expand(
            f"{TAXONOMY_ASSEMBLY_DIR}/{{sample}}/{wc.marker}/alignment/{{sample}}.idxstats.tsv",
            sample=SAMPLES,
        )
    output:
        matrix=f"{TAXONOMY_ASSEMBLY_DIR}/summary/{{marker}}/matrices/contig_counts.tsv"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.matrix}")"
        awk '
            BEGIN {{ FS = OFS = "\t" }}

            FNR == 1 {{
                file++
                row = 0
                sample = FILENAME
                sub(/^.*\//, "", sample)
                sub(/\.idxstats\.tsv$/, "", sample)
                samples[file] = sample
            }}

            $1 == "*" {{ next }}

            {{
                row++
                if (file == 1) {{
                    contigs[row] = $1
                    n_contigs = row
                }} else if (row > n_contigs || contigs[row] != $1) {{
                    print "Contig order or IDs differ in " FILENAME > "/dev/stderr"
                    failed = 1
                    exit 1
                }}
                counts[row, file] = $3
                rows[file] = row
            }}

            END {{
                if (failed) exit 1
                if (file != ARGC - 1 || n_contigs == 0) {{
                    print "Missing or empty idxstats input" > "/dev/stderr"
                    exit 1
                }}
                for (i = 1; i <= file; i++) {{
                    if (rows[i] != n_contigs) {{
                        print "Contig count differs for " samples[i] > "/dev/stderr"
                        exit 1
                    }}
                }}

                printf "%s", "ContigID"
                for (i = 1; i <= file; i++) printf "\t%s", samples[i]
                printf "\n"

                for (row = 1; row <= n_contigs; row++) {{
                    printf "%s", contigs[row]
                    for (i = 1; i <= file; i++) printf "\t%s", counts[row, i]
                    printf "\n"
                }}
            }}
        ' {input.stats:q} > "{output.matrix}"
        """