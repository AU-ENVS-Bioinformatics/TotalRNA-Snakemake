from pathlib import Path
from snakemake.io import expand


def get_hmmscan_outputs(wildcards):

    checkpoint_output = checkpoints.hmmscan_split.get()
    chunk_dir = Path(checkpoint_output.output.chunks)

    chunk_files = sorted(chunk_dir.glob("*.pep"))

    print(f"[DEBUG] found {len(chunk_files)} chunks")

    if not chunk_files:
        raise ValueError(
            f"No chunk files found in {chunk_dir}"
        )

    return expand(
        f"{RESULTS_DIR}/nonrRNA/hmmscan/results/{{chunk}}.pfam.domtblout",
        chunk=[f.stem for f in chunk_files]
    )

checkpoint hmmscan_split:
    conda:
        "../envs/seqkit.yaml"
    input:
        pep=f"{RESULTS_DIR}/nonrRNA/predicted/longest_orfs.pep"
    output:
        chunks=directory(
            f"{RESULTS_DIR}/nonrRNA/hmmscan/chunks"
        )
    params:
        chunk_size=config["nonrRNA"]["hmmscan"]["chunk_size"]
    shell:
        r"""
        rm -rf {output.chunks}

        mkdir -p {output.chunks}

        seqkit split2 \
            -t protein \
            -s {params.chunk_size} \
            -O {output.chunks} \
            {input.pep}
        """


rule merge_hmmscan:
    input:
        get_hmmscan_outputs
    output:
        pfam=f"{RESULTS_DIR}/nonrRNA/hmmscan/pfam.domtblout"
    run:
        Path(output.pfam).parent.mkdir(
            parents=True,
            exist_ok=True
        )
        first = True
        with open(output.pfam, "w") as out:
            for infile in input:
                with open(infile) as fin:
                    for line in fin:
                        if line.startswith("#"):
                            if first:
                                out.write(line)
                        else:
                            out.write(line)
                first = False