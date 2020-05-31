library(data.table)

args <- commandArgs(trailingOnly = TRUE)
message("Input:")
print(args)

if(length(args) == 1) {
  config <- args[1]
} else {
  stop("Incorrect arguments")
}
rm(args)

source(config)
source("gene_type.R")

## Mean exonic coverage for each sample
cov.path <- paste0(covDIR, "/merge/exon/")
mean.cov <- function(genes) {
    cov <- data.table()
    for (g in genes) {
        cov <- rbind(cov, fread(paste0(cov.path, g, "_exon.cov"), sep=" ", header=TRUE))
        # The input file for each gene contains a table with coverage of each exonic base for each sample.
        # Structure of this file (space-separated):
        # loc           Sample1 Sample2 Sample3 Sample4 ...
        # start_exon1   174     168     268     240     ...
        # start_exon1+1 175     175     265     242     ...
        # ...
        # end_exon1     184     177     271     247     ...
        # start_exon2   178     180     270     240     ...
        # ...
    }
    apply(cov[, -1], 2, mean)
}

cat("Calculating mean coverage for this batch\n", file=stdout())
sample.mean.cov <- data.frame(t(rbind(mean.cov(auto.genes), mean.cov(X.genes))))
colnames(sample.mean.cov) <- c("auto.genes", "sex.genes")
print(sample.mean.cov)
save(sample.mean.cov, file=paste0(CNVDIR, "/", BATCH, "-mean-cov.rda"))
