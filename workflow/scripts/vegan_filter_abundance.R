run <- function(infile, outfile, minimum) {
    # Fail if minimum is not a number
    if (!is.numeric(minimum)) {
        stop("Minimum must be a number")
    }
    # Fail if vegan is not installed
    if (!requireNamespace("vegan", quietly = TRUE)) {
        stop("vegan package is not installed")
    }
    data <- read.table(infile,sep="\t",header=TRUE,row.names=1)
    lowest_sample_reads_sum <- data |> colSums() |> min()
    data_stand <- data |> t() |> vegan::decostand(method="total") |> t()
    included_idx <- rowMeans(data_stand) >= minimum/lowest_sample_reads_sum
    included_contigs <- row.names(data_stand)[included_idx]
    # Data rownames to column ContigID as first column using base R
    data$ContigID <- row.names(data)
    data <- data[,c(ncol(data),1:(ncol(data)-1))]
    write.table(data[included_contigs,], outfile, sep="\t", quote=FALSE, row.names=FALSE)

}


run(snakemake@input[[1]], snakemake@output[[1]], snakemake@params[["minimum"]])