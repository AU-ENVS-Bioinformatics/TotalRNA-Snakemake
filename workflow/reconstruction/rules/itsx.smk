rule ITSx:
    conda:
        "../envs/itsx.yaml"
    message:
        "[ITSx] Extracting ITS regions from the reconstructed ITS candidates for {wildcards.sample}"
    input:
        transcripts=f"{ASSEMBLY_DIR}/cross_assembly/{{sample}}/ITS_candidates/rnaspades/transcripts.fasta"
    output:
        its1=f"{ASSEMBLY_DIR}/cross_assembly/{{sample}}/ITS_candidates/ITSx/{{sample}}_ITSx.ITS1.fasta",
        five_8s=f"{ASSEMBLY_DIR}/cross_assembly/{{sample}}/ITS_candidates/ITSx/{{sample}}_ITSx.5_8S.fasta",
        its2=f"{ASSEMBLY_DIR}/cross_assembly/{{sample}}/ITS_candidates/ITSx/{{sample}}_ITSx.ITS2.fasta",
        full=f"{ASSEMBLY_DIR}/cross_assembly/{{sample}}/ITS_candidates/ITSx/{{sample}}_ITSx.full.fasta",
        its1_partial=f"{ASSEMBLY_DIR}/cross_assembly/{{sample}}/ITS_candidates/ITSx/{{sample}}_ITSx.ITS1.full_and_partial.fasta",
        five_8s_partial=f"{ASSEMBLY_DIR}/cross_assembly/{{sample}}/ITS_candidates/ITSx/{{sample}}_ITSx.5_8S.full_and_partial.fasta",
        its2_partial=f"{ASSEMBLY_DIR}/cross_assembly/{{sample}}/ITS_candidates/ITSx/{{sample}}_ITSx.ITS2.full_and_partial.fasta",
        full_partial=f"{ASSEMBLY_DIR}/cross_assembly/{{sample}}/ITS_candidates/ITSx/{{sample}}_ITSx.full_and_partial.fasta",
        summary=f"{ASSEMBLY_DIR}/cross_assembly/{{sample}}/ITS_candidates/ITSx/{{sample}}_ITSx.summary.txt"
    log:
        stdout=f"{ASSEMBLY_DIR}/cross_assembly/{{sample}}/ITS_candidates/ITSx/{{sample}}_ITSx.log"
    benchmark:
        f"{ASSEMBLY_DIR}/cross_assembly/{{sample}}/ITS_candidates/ITSx/{{sample}}_ITSx_benchmark.txt"
    params:
        options=config["reconstruction"]["itsx"].get("options", ""),
        species=config["reconstruction"]["itsx"].get("species", "Fungi"),
        regions=config["reconstruction"]["itsx"].get("regions", "ITS1,5.8S,ITS2"),
        output_prefix=f"{ASSEMBLY_DIR}/cross_assembly/{{sample}}/ITS_candidates/ITSx/{{sample}}_ITSx"
    threads:
        config["reconstruction"]["itsx"].get("threads", 12)
    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {params.output_prefix})

        ITSx \
            -i {input.transcripts} \
            -o {params.output_prefix} \
            -t {params.species} \
            --cpu {threads} \
            --save_regions {params.regions} \
            {params.options} \
            > {log.stdout} 2>&1
        
        touch {output.its1} {output.five_8s} {output.its2} {output.full} {output.its1_partial} {output.five_8s_partial} {output.its2_partial} {output.full_partial}
        """