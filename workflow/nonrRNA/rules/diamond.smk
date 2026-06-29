nonrrna_module = config["nonrRNA"]["nonrrna_module"]

if nonrrna_module in ["coassembly", "both", "all"]:
    if config["nonrRNA"]["assembly_predictor"] == "transdecoder":
        rule diamond_blastp:
            conda:
                "../envs/diamond.yaml"
            message:
                "[DIAMOND] protein alignment (UniRef90)"
            input:
                pep=f"{RESULTS_DIR}/nonrRNA/predicted/longest_orfs.pep"
            output:
                diamond_uniref90=f"{RESULTS_DIR}/nonrRNA/diamond/transdecoder.blastp.outfmt6",
            log:
                stdout=f"{RESULTS_DIR}/nonrRNA/diamond/logs/diamond.log"
            benchmark:
                f"{RESULTS_DIR}/nonrRNA/diamond/benchmarks/diamond.txt"
            params:
                db=config["databases"]["diamond_proteindb"],
                options=config["nonrRNA"]["diamond_blastp"]["options"]
            threads:
                config["nonrRNA"]["diamond_blastp"].get("threads", 14)
            shell:
                r"""
                set -euo pipefail
                mkdir -p $(dirname {output.diamond_uniref90})

                diamond blastp \
                    -q {input.pep} \
                    -d {params.db} \
                    -o {output.diamond_uniref90} \
                    --threads {threads} \
                    {params.options} \
                    > {log.stdout} 2>&1
                """

#diamond blastp -q longest_orfs.pep 
#--db /data_2/Databases/uniref90_db/uniref90.dmnd 
#--out transdecoder.blastp.outfmt6 
#--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore 
#--evalue 1e-5 --max-target-seqs 1 --sensitive --threads 16