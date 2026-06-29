nonrrna_module = config["nonrRNA"]["nonrrna_module"]

if nonrrna_module in ["coassembly", "both", "all"]:
    if config["nonrRNA"]["assembly_method"] in ["spades","rnaspades"]:
        rule rnaSpades:
            conda:
                "../envs/spades.yaml"
            message:
                "[rnaSpades] co-assembly across all samples"
            input:
                concatenated_fastq_r1 = f"{RESULTS_DIR}/nonrRNA/concatenated/nonrRNA_1.fastq.gz",
                concatenated_fastq_r2 = f"{RESULTS_DIR}/nonrRNA/concatenated/nonrRNA_2.fastq.gz"
            output:
                assembly=f"{RESULTS_DIR}/nonrRNA/coassembly/transcripts.fasta"
            log:
                stdout=f"{RESULTS_DIR}/nonrRNA/coassembly/logs/rnaspades.log"
            benchmark:
                f"{RESULTS_DIR}/nonrRNA/coassembly/benchmarks/rnaspades.txt"
            threads:
                config["nonrRNA"]["rnaspades"].get("threads", 8)
            params:
                options=config["nonrRNA"]["rnaspades"]["options"],
                outdir=f"{RESULTS_DIR}/nonrRNA/coassembly/"
            shell:
                r"""
                set -euo pipefail

                mkdir -p {params.outdir}
                mkdir -p $(dirname {log})

                rnaspades.py \
                    -1 {input.concatenated_fastq_r1} \
                    -2 {input.concatenated_fastq_r2} \
                    -t {threads} \
                    -o {params.outdir} \
                    {params.options}
                    > {log.stdout} 2>&1
                """
    # in the future additional assemblers, like trinity can be reintroduced or metaspades or metahit
