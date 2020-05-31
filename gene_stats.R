library(GenomicRanges)

args <- commandArgs(trailingOnly = TRUE)
message("Input:")
print(args)

if(length(args) == 2) {
  config <- args[1]
  gene <- args[2]
} else {
  stop("Incorrect arguments")
}
rm(args)

source(config)

all.genes <- read.table(paste0(PreProcessDIR, "/", PROJECT, "-genes.bed"), sep="\t", header=FALSE, col.names=c("chrom", "start", "end", "gene"))
chr <- paste0("chr", all.genes[all.genes$gene == gene, "chrom"][1])

gene.raw <- read.table(paste0(covDIR, "/merge/gene/", gene, ".cov"), sep=" ", fill=TRUE, header=TRUE)
# The input file for each gene contains a table with coverage of all gene's bases (including intronic) with the flanks for each sample.
# Structure of this file (space-separated):
# loc                Sample1 Sample2 Sample3 Sample4 ...
# start_gene-flank   112     94      171     140     ...
# start_gene-flank+1 114     96      175     145     ...
# ...
# end_gene+flank     119     98      185     154     ...
gene.raw[is.na(gene.raw)] <- 0

gene.raw.noloc <- data.frame(gene.raw[-1])
gene.raw.loc <- data.frame(gene.raw[1])
gene.cov.median <- apply(gene.raw.noloc, 1, median, na.rm = TRUE)
gene.cov.lower <- apply(gene.raw.noloc, 1, quantile, probs = c(0.05),  na.rm = TRUE)
gene.cov.upper <- apply(gene.raw.noloc, 1, quantile, probs = c(0.95),  na.rm = TRUE)

gene.cov.summary <- data.frame(gene.raw.loc, gene.cov.lower, gene.cov.median, gene.cov.upper)
colnames(gene.cov.summary) <- c("loc", "lower", "median", "upper")

df.gr.gene.cov.summary <- cbind(data.frame(chr=chr, start=gene.cov.summary$loc, end=gene.cov.summary$loc,
                                           strand="*"), gene.cov.summary[,2:4])
save(df.gr.gene.cov.summary, file = paste0(covDIR, "/geneplot/", gene, ".rda"))
