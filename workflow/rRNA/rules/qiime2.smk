#https://docs.qiime2.org/2024.10/tutorials/importing/ -> “Fastq manifest” formats

rule qiime_manifest:
    conda: "../envs/python.yaml"
    message:
        "[QIIME2] generate manifest for stage {wildcards.stage}"
    input:
        aligned_fwd=expand(
            f"{RESULTS_DIR}/qc/{{sample}}/sortmerna/{{stage}}/{{sample}}_{{stage}}.aligned_fwd.fq.gz",
            sample=SAMPLES,
            stage=lambda wc: wc.stage
        ),
        aligned_rev=expand(
            f"{RESULTS_DIR}/qc/{{sample}}/sortmerna/{{stage}}/{{sample}}_{{stage}}.aligned_rev.fq.gz",
            sample=SAMPLES,
            stage=lambda wc: wc.stage
        )
    log:
        stdout=f"{RESULTS_DIR}/rRNA/logs/qiime_{{stage}}_manifest.log"
    benchmark:
        f"{RESULTS_DIR}/rRNA/benchmarks/qiime_{{stage}}_manifest.txt"
    output:
        manifest=f"{RESULTS_DIR}/rRNA/qiime/{{stage}}/qiime_{{stage}}_manifest.tsv"
    run:
        import os

        os.makedirs(os.path.dirname(output.manifest), exist_ok=True)

        with open(output.manifest, "w") as fh:
            fh.write(
                "sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n"
            )

            for sample, fwd, rev in zip(
                SAMPLES,
                input.aligned_fwd,
                input.aligned_rev
            ):
                fh.write(
                    f"{sample}\t{os.path.abspath(fwd)}\t{os.path.abspath(rev)}\n"
                )


#https://docs.qiime2.org/2024.10/tutorials/importing/
rule qiime_import:
    conda: "../envs/qiime2.yaml"
    message:
        "[QIIME2] import reads for across samples for stage {wildcards.stage}"
    input:
        manifest=f"{RESULTS_DIR}/rRNA/qiime/{{stage}}/qiime_{{stage}}_manifest.tsv"
    output:
        qza=f"{RESULTS_DIR}/rRNA/qiime/{{stage}}/qiime_{{stage}}_demux.qza"
    log:
        stdout = f"{RESULTS_DIR}/rRNA/logs/qiime_{{stage}}_import.log"
    benchmark:
        f"{RESULTS_DIR}/rRNA/benchmarks/qiime_{{stage}}_import.txt"
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.qza})

        qiime tools import \
        --type 'SampleData[PairedEndSequencesWithQuality]' \
        --input-format PairedEndFastqManifestPhred33 \
        --input-path {input.manifest} \
        --output-path {output.qza}
        """
