## Gene types
# All input files are lists of gene names, one gene per line.
# We always have autosomal and non-PAR X-chromosome genes.
# We don't have diagnostic genes on the Y chromosome.
# We may have genes in the PAR of the X chromosome.
# The genes on the PAR of the X chromosome are treated the same way as the autosomal genes.

auto.genes <- read.table(paste0(PreProcessDIR, "/", PROJECT, "-AUTOgenes"), stringsAsFactors=FALSE)$V1

options(show.error.messages = FALSE)
par.genes <- try(read.table(paste0(PreProcessDIR, "/", PROJECT, "-PARgenes"), stringsAsFactors=FALSE))
if(! inherits(par.genes, "try-error")) {
    auto.genes <- c(auto.genes, par.genes$V1)
}
options(show.error.messages = TRUE)

X.genes <- read.table(paste0(PreProcessDIR, "/", PROJECT, "-Xgenes"), stringsAsFactors=FALSE)$V1  # those are genes in the non-PAR of the X
