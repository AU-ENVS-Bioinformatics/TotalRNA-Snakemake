################################################################################
# HMMscan PFAM EVIDENCE FOR TRANSDECODER
################################################################################

FUNCTIONAL_MODULE = config["functional_profiling"]["module"]
ORF_PREDICTOR = config["functional_profiling"]["ORF_predictor"]

if FUNCTIONAL_MODULE in ["assembly", "both", "all"] and ORF_PREDICTOR == "transdecoder":

    rule hmmscan_chunk:
        conda:
            "../envs/hmmer.yaml"
        message:
            "[HMMscan] Search peptide chunk {wildcards.chunk} against Pfam"
        input:
            pep=f"{FUNCTION_DIR}/hmmscan/chunks/{{chunk}}.pep"
        output:
            domtblout=f"{FUNCTION_DIR}/hmmscan/results/{{chunk}}.pfam.domtblout",
            report=f"{FUNCTION_DIR}/hmmscan/results/{{chunk}}.hmmscan.out"
        log:
            stdout=f"{FUNCTION_DIR}/hmmscan/logs/{{chunk}}.hmmscan.log",
            time=f"{FUNCTION_DIR}/hmmscan/benchmarks/{{chunk}}.hmmscan_time.txt"
        benchmark:
            f"{FUNCTION_DIR}/hmmscan/benchmarks/{{chunk}}.hmmscan.txt"
        params:
            db=config["databases"]["Pfam_db"],
            options=config["functional_profiling"]["hmmscan"]["options"]
        threads:
            config["functional_profiling"]["hmmscan"].get("threads", 2)
        shell:
            r"""
            set -euo pipefail

            mkdir -p $(dirname {output.domtblout})

            /usr/bin/time -v -o {log.time} \
                hmmscan \
                    --cpu {threads} \
                    --domtblout {output.domtblout} \
                    {params.options} \
                    -o {output.report} \
                    {params.db} \
                    {input.pep} \
                    > {log.stdout} 2>&1
            """