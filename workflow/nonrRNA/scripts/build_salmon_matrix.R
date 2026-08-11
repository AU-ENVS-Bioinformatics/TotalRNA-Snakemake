#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
    stop(
        paste(
            "Usage:",
            "Rscript build_salmon_matrix.R",
            "<column_name>",
            "<output_file>",
            "<quant.sf> [quant.sf ...]"
        )
    )
}

value_col <- args[1]
output_file <- args[2]
quant_files <- args[-c(1, 2)]

outdir <- dirname(output_file)

if (!dir.exists(outdir) && outdir != ".") {
    dir.create(
        outdir,
        recursive = TRUE
    )
}

message(
    "Found ",
    length(quant_files),
    " quant.sf files"
)

sample_tables <- lapply(
    quant_files,
    function(file) {

        sample_name <- basename(
            dirname(file)
        )

        message(
            "Reading ",
            sample_name
        )

        dt <- fread(file)

        if (!value_col %in% colnames(dt)) {
            stop(
                "Column '",
                value_col,
                "' not found in ",
                file
            )
        }

        dt <- dt[
            ,
            .(
                Name,
                value = get(value_col)
            )
        ]

        setnames(
            dt,
            "value",
            sample_name
        )

        dt
    }
)

count_matrix <- Reduce(
    function(x, y) {
        merge(
            x,
            y,
            by = "Name",
            all = TRUE
        )
    },
    sample_tables
)

count_matrix[is.na(count_matrix)] <- 0

fwrite(
    count_matrix,
    file = output_file,
    sep = "\t"
)

message(
    "Matrix dimensions: ",
    nrow(count_matrix),
    " features x ",
    ncol(count_matrix) - 1,
    " samples"
)

message(
    "Output written to ",
    normalizePath(output_file)
)