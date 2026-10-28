################################################################################
# IDENTIFY DOMAINS WITH NONEMPTY FINAL SSU FASTAS
################################################################################

checkpoint check_phylogeny_domains:
    input:
        fasta=expand(f"{ASSEMBLY_DIR}/cross_assembly/SSU/domains/{{domain}}/final/transcript.fasta", domain=PHYLOGENY_DOMAINS)
    output:
        status=f"{PHYLOGENY_DIR}/domain_status.tsv"
    run:
        Path(output.status).parent.mkdir(parents=True, exist_ok=True)

        with open(output.status, "w") as report:
            report.write("domain\tstatus\n")

            for domain, fasta in zip(PHYLOGENY_DOMAINS, input.fasta):
                with open(fasta) as handle:
                    present = any(line.startswith(">") for line in handle)

                status = "present" if present else "empty"
                report.write(f"{domain}\t{status}\n")