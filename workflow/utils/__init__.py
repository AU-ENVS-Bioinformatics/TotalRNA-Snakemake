"""
Module-specific helper functions for the TotalRNA Snakemake workflow.

This package separates reusable Python logic from the main Snakefile.
The Snakefile coordinates the workflow, while these modules provide
functions for selecting inputs, validating settings, and constructing
paths.

Planned helper modules:
    qc:
        Raw-read lookup and sample-specific QC directories.
    rna_classification:
        RNA input selection, SortMeRNA stages, and reference lookup.
    reconstruction:
        Reconstruction settings and Barrnap status-file parsing.
    taxonomy:
        Taxonomy input selection, database lookup, and tool settings.
    functional_profiling:
        Helpers for read-based and assembly-based functional analysis.
    phylogeny:
        Reference selection and BLAST status-file parsing.
    habitat_refinement:
        Habitat-specific inputs, databases, and configuration validation.
    viral_sequence_analysis:
        Helpers for viral identification and downstream analysis.

Only modules needed by the workflow need to be implemented.

The targets package defines the outputs required to complete each
workflow module. This utils package provides the supporting logic used
to obtain those outputs.

Helper functions receive configuration, sample information, and paths
explicitly rather than relying on variables defined in the Snakefile.
This makes their dependencies clear and allows them to be reused and
checked independently.

This __init__.py file marks utils as a regular Python package, allowing
explicit imports such as `from utils import qc`. It does not execute
workflow steps or load configuration. Modern Python can support imports
without this file, but keeping it makes the package structure explicit.
"""