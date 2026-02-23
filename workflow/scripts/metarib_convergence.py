#!/usr/bin/env python3
"""
Check if MetaRib iteration should continue or stop.

This script evaluates convergence criteria for the iterative assembly process:
1. Minimum unmapped reads threshold
2. Estimated alignment rate (predict EMIRGE failure)
3. Read depletion rate
4. Contig assembly convergence (absolute and relative)

Used as Snakemake script directive (receives snakemake object).
Writes decision to output file (True = STOP, False = CONTINUE).
"""
from pathlib import Path


def read_contig_read_count(report_path):
    """Read contig and unmapped read counts from report file.
    
    Returns:
        tuple: (contig_count, unmapped_reads)
    """
    path = Path(report_path)
    if not path.exists():
        return 0, 0

    with open(path) as f:
        contig_count = int(f.readline().strip() or 0)
        unmapped_reads = int(f.readline().strip() or 0)
        return contig_count, unmapped_reads


def check_convergence(
    current_report,
    previous_report,
    min_reads_threshold,
    convergence_threshold,
    log_file,
    iteration=None,
):
    """
    Check if iteration should continue based on convergence criteria.
    
    Returns:
        bool: True if should STOP, False if should CONTINUE
    """
    log_lines = []
    
    if iteration is not None:
        log_lines.append(f"Checking iteration {iteration}")

    # Read metrics from reports
    contig_prev, unmapped_prev = read_contig_read_count(previous_report)
    contig_curr, unmapped_curr = read_contig_read_count(current_report)
    
    log_lines.extend([
        f"Previous iteration: {contig_prev} contigs, {unmapped_prev} unmapped reads",
        f"Current iteration: {contig_curr} contigs, {unmapped_curr} unmapped reads",
    ])

    # Calculate metrics
    depletion_rate = (unmapped_prev - unmapped_curr) / unmapped_prev if unmapped_prev > 0 else 0
    estimated_alignment_rate = depletion_rate * 0.5  # heuristic
    contig_diff = contig_curr - contig_prev
    
    log_lines.extend([
        f"Read depletion rate: {depletion_rate:.2%}",
        f"Estimated alignment rate for next iteration: {estimated_alignment_rate:.2%}",
        f"Absolute contig change: {contig_diff:+d}",
    ])

    # Determine convergence with early returns for clarity
    def stop(reason):
        """Helper to indicate stopping with reason."""
        log_lines.extend(["", "DECISION: STOP", f"REASON: {reason}"])
        if log_file:
            Path(log_file).write_text("\n".join(log_lines) + "\n")
        return True

    def continue_iter(reason):
        """Helper to indicate continuing with reason."""
        log_lines.extend(["", "DECISION: CONTINUE", f"REASON: {reason}"])
        if log_file:
            Path(log_file).write_text("\n".join(log_lines) + "\n")
        return False

    # Check stopping criteria
    if unmapped_curr <= min_reads_threshold:
        return stop(
            f"Unmapped reads below threshold ({unmapped_curr} <= {min_reads_threshold})"
        )
    
    if contig_curr == 0:
        return stop("No contigs assembled in current iteration")
    
    if estimated_alignment_rate < 0.05:
        return stop(
            f"Predicted alignment rate too low ({estimated_alignment_rate:.2%} < 5%). "
            "EMIRGE likely to fail due to insufficient coverage."
        )
    
    if depletion_rate < 0.10 and unmapped_prev > 0:
        return stop(
            f"Read depletion stalled ({depletion_rate:.2%} < 10%). "
            "Most rRNA sequences likely assembled."
        )

    # Check contig convergence
    if contig_prev > 0:
        if abs(contig_diff) <= 1:
            return stop(f"Contig count converged (Δ={contig_diff:+d}). Assembly stable.")
        
        if contig_diff > 0:
            pct_change = contig_diff / contig_prev
            log_lines.append(f"Relative contig change: {pct_change:.2%}")
            
            if pct_change <= convergence_threshold:
                return stop(
                    f"Contig growth below threshold ({pct_change:.2%} <= {convergence_threshold:.2%}). "
                    "Convergence reached."
                )
            else:
                return continue_iter(
                    f"Significant contig growth detected ({pct_change:.2%} > {convergence_threshold:.2%})"
                )
        else:
            return stop(f"Contig count decreased ({contig_diff:+d}). Converged.")
    
    # First iteration with contigs
    return continue_iter("First iteration with contigs")


# Main execution when called by Snakemake
should_stop = check_convergence(
    current_report=snakemake.input.cur_report,
    previous_report=snakemake.input.prev_report,
    min_reads_threshold=snakemake.params.min_reads_threshold,
    convergence_threshold=snakemake.params.convergence_threshold,
    log_file=snakemake.log[0],
    iteration=snakemake.wildcards.iter,
)

# Write decision to output file
Path(snakemake.output[0]).write_text("True\n" if should_stop else "False\n")

