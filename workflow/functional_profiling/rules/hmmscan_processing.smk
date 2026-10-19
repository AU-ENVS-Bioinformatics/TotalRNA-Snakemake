from pathlib import Path

################################################################################
# SPLIT TRANSDECODER CANDIDATE PROTEINS
################################################################################

checkpoint hmmscan_split:
    conda:
        "../envs/seqkit.yaml"
    message:
        "[HMMscan] Split candidate proteins into chunks"
    input:
        pep=f"{FUNCTION_DIR}/ORF_prediction/transdecoder/longorfs/longest_orfs.pep"
    output:
        chunks=directory(f"{FUNCTION_DIR}/hmmscan/chunks")
    log:
        stdout=f"{FUNCTION_DIR}/hmmscan/logs/split.log"
    benchmark:
        f"{FUNCTION_DIR}/hmmscan/benchmarks/split.txt"
    params:
        chunk_size=config["functional_profiling"]["hmmscan"]["chunk_size"]
    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {output.chunks}) 

        seqkit split2 \
            -t protein \
            -s {params.chunk_size} \
            -O {output.chunks} \
            {input.pep} \
            > {log.stdout} 2>&1
        """

################################################################################
# MERGE HMMscan RESULTS
################################################################################

rule merge_hmmscan:
    message:
        "[HMMscan] Merge Pfam domain results"
    input:
        domtblout=lambda wc: get_hmmscan_outputs(
            checkpoints=checkpoints,
            results_dir=RESULTS_DIR,
        )
    output:
        pfam=f"{FUNCTION_DIR}/hmmscan/merged/pfam.domtblout"
    run:
        Path(output.pfam).parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        first_file = True

        with open(output.pfam, "w") as outfile:
            for infile in input.domtblout:
                with open(infile) as source:
                    for line in source:
                        if line.startswith("#"):
                            if first_file:
                                outfile.write(line)
                        else:
                            outfile.write(line)

                first_file = False