from pathlib import Path
from snakemake.io import expand


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
# DISCOVER HMMscan CHUNK OUTPUTS
################################################################################

def get_hmmscan_outputs(wildcards):

    checkpoint_output = checkpoints.hmmscan_split.get()
    chunk_dir = Path(checkpoint_output.output.chunks)

    chunk_files = sorted(chunk_dir.glob("*.pep"))

    if not chunk_files:
        raise ValueError(
            f"No peptide chunks found in {chunk_dir}"
        )

    return expand(
        f"{FUNCTION_DIR}/hmmscan/results/{{chunk}}.pfam.domtblout",
        chunk=[chunk.stem for chunk in chunk_files],
    )


################################################################################
# MERGE HMMscan RESULTS
################################################################################

rule merge_hmmscan:
    message:
        "[HMMscan] Merge Pfam domain results"
    input:
        domtblout=get_hmmscan_outputs
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