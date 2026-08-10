nonrrna_module = config["nonrRNA"]["nonrrna_module"]

if nonrrna_module in ["coassembly", "both", "all"]:
    if config["nonrRNA"]["assembly_predictor"] == "transdecoder":

        rule hmmscan_chunk:
            conda:
                "../envs/hmmer.yaml"
            input:
                pep=f"{RESULTS_DIR}/nonrRNA/hmmscan/chunks/{{chunk}}.pep"
            output:
                domtblout=f"{RESULTS_DIR}/nonrRNA/hmmscan/results/{{chunk}}.pfam.domtblout",
                direct_output=f"{RESULTS_DIR}/nonrRNA/hmmscan/results/{{chunk}}.hmmscan.out"
            log:
                stdout=f"{RESULTS_DIR}/nonrRNA/hmmscan/logs/{{chunk}}.hmmscan.log",
                benchmark_file=f"{RESULTS_DIR}/nonrRNA/hmmscan/benchmarks/{{chunk}}.hmmscan_time.txt",
            benchmark:
                f"{RESULTS_DIR}/nonrRNA/hmmscan/benchmarks/{{chunk}}.hmmscan.txt"
            params:
                db=config["databases"]["Pfam_db"],
                options=config["nonrRNA"]["hmmscan"]["options"]
            threads:
                config["nonrRNA"]["hmmscan"].get("threads", 2)
            shell:
                r"""
                mkdir -p $(dirname {output.domtblout})
                mkdir -p $(dirname {log.stdout})

                /usr/bin/time -v -o {log.benchmark_file} \
                    hmmscan \
                        --cpu {threads} \
                        --domtblout {output.domtblout} \
                        {params.options} \
                        -o {output.direct_output} \
                        {params.db} \
                        {input.pep} \
                        > {log.stdout} 2>&1
                """