if config["nonrRNA"]["module"] == "coassembly" and config["nonrRNA"]["predictor"] == "transdecoder":

    rule hmmscan:
        conda:
            "../envs/hmmer.yaml"
        message:
            "[HMMER] Pfam domain scan"
        input:
            pep=f"{RESULTS_DIR}/nonrRNA/predicted/longest_orfs.pep"
        output:
            pfam=f"{RESULTS_DIR}/nonrRNA/hmmscan/pfam.domtblout",
            logfile=f"{RESULTS_DIR}/nonrRNA/hmmscan/stdout.txt"
        log:
            stdout = f"{RESULTS_DIR}/nonrRNA/hmmscan/logs/hmmscan.log"
        benchmark:
            f"{RESULTS_DIR}/nonrRNA/hmmscan/benchmarks/hmmscan.txt"
        params:
            db=config["databases"]["Pfam_db"]
        threads:
            config["nonrRNA"]["hmm_scan"].get("threads", 16)
        shell:
            r"""
            set -euo pipefail
            mkdir -p $(dirname {output.pfam})

            hmmscan \
                --cpu {threads} \
                --domtblout {output.pfam} \
                -o {output.logfile} \
                {params.db} \
                {input.pep} \
                > {log.stdout} 2>&1
            """