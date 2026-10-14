################################################################################
# CHECK WHICH SSU CATALOGUES CONTAIN SEQUENCES
################################################################################

checkpoint kingdom_check:
    message:
        "[SSU] Checking kingdom catalogues before phylogenetic analysis"
    input:
        catalogues=expand(f"{ASSEMBLY_DIR}/SSU_candidates/postprocessing/{{kingdom}}/SSU_{{kingdom}}_catalogue.fasta", kingdom=BARRNAP_KINGDOMS)
    output:
        report=f"{ASSEMBLY_DIR}/SSU_candidates/postprocessing/kingdom_check.tsv"
    run:
        from pathlib import Path

        Path(output.report).parent.mkdir(parents=True, exist_ok=True)

        with open(output.report, "w") as report:
            report.write("kingdom\tsequence_count\tstatus\n")

            for kingdom, fasta_path in zip(BARRNAP_KINGDOMS, input.catalogues):
                sequence_count = 0
                current_id = None
                current_length = 0

                with open(fasta_path) as handle:
                    for line in handle:
                        if line.startswith(">"):
                            if current_id is not None and current_length == 0:
                                raise ValueError(f"Empty FASTA record {current_id} in {fasta_path}")

                            header = line[1:].strip()

                            if not header:
                                raise ValueError(f"Empty FASTA header in {fasta_path}")

                            current_id = header.split()[0]
                            current_length = 0
                            sequence_count += 1

                        elif line.strip():
                            if current_id is None:
                                raise ValueError(f"Sequence before FASTA header in {fasta_path}")

                            current_length += len("".join(line.split()))

                if current_id is not None and current_length == 0:
                    raise ValueError(f"Empty FASTA record {current_id} in {fasta_path}")

                status = "available" if sequence_count else "empty"
                report.write(f"{kingdom}\t{sequence_count}\t{status}\n")
                print(f"[SSU] {kingdom}: {sequence_count} sequences; {status}")