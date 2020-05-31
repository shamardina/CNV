library(data.table)
library(ggplot2)
library(Biobase)
library(biomaRt)
library(Gviz)
library(GenomicRanges)

args <- commandArgs(trailingOnly = TRUE)
message("Input:")
print(args)

if(length(args) == 3) {
  config <- args[1]
  sample <- args[2]
  gene <- args[3]
} else {
  stop("Incorrect arguments")
}
rm(args)

source(config)

## ------------------------------------------------------------------------
## Tracks
ref <- "hg19"
df.header <- c("chrom", "start", "end", "gene")

## Design gene regions (based on the chosen transcript):
# 2000bp before the gene's start and 100bp after the gene's end (taking strand into account)
design <- fread(paste0(PreProcessDIR, "/", PROJECT, "-genes-design.bed"), sep="\t", data.table=FALSE, header=FALSE, col.names = df.header)
window.start <- design[design$gene == gene, "start"] + 1 - 1000
window.end <- design[design$gene == gene, "end"] + 1000
chr <- paste0("chr", design[design$gene == gene, "chrom"][1])

## The gene's panel target
targets <- fread(paste0(PreProcessDIR, "/", PROJECT, "-target-flank0.bed"), sep="\t", data.table=FALSE, header=FALSE, col.names = df.header[1:3])
targets$start <- targets$start + 1

## Problematic regions:
# regions with coverage <20X in at least 5% of the samples
# (see Simeoni I, Shamardina O, Deevi SV, et al. GRID – Genomics of Rare Immune Disorders...)
pr5p <- fread(paste0(covDIR, "/", PROJECT, "-problematic.5p.interval"), sep="\t", data.table=FALSE, header=FALSE, col.names = df.header)
pr5p <- pr5p[pr5p$chrom == design[design$gene == gene, "chrom"][1] & pr5p$start >= window.start & pr5p$end <= window.end, ]

## Ensembl transcript ID
transcripts <- fread(paste0(PreProcessDIR, "/", PROJECT, "-transcript-gene.info"), sep="\t", data.table=FALSE, header=FALSE, col.names = c("transcriptid", "symbol"))
transcript.id <- transcripts[transcripts$symbol == gene,]$transcriptid

## CNV calling regions
CNVregions <- fread(paste0(CNVDIR, "/cnv-regions.bed"), sep="\t", data.table=FALSE, header=FALSE, col.names=df.header)
CNVregions$start <- CNVregions$start + 1
CNVregions <- CNVregions[CNVregions$chrom == design[design$gene == gene, "chrom"][1] & CNVregions$start >= window.start & CNVregions$end <= window.end, ]
print(CNVregions)

## CNV problematic regions:
# the regions with median relative normalized read counts < 0.2
# (see Simeoni I, Shamardina O, Deevi SV, et al. GRID – Genomics of Rare Immune Disorders...)
# The header of the file is:
# space   start   end     width   Min.    Median  Max.    Batch
CNVprregions <- fread(paste0(CNVDIR, "/out/poorly_covered_regions.tsv"), sep="\t", data.table=FALSE, header = TRUE)
CNVprregions <- CNVprregions[CNVprregions$space == design[design$gene == gene, "chrom"][1] & CNVprregions$start >= window.start & CNVprregions$end <= window.end, ]

## Sample's CNV in this gene:
# This table contains merged and annotated results of CNV calling.
# The selected fields from the table are:
# sample  type  id  BF  genes  promoters
# 'type' is 'deletion' or 'duplication'
# 'BF' is Bayes factor reported by the CNV caller (ExomeDepth R-package version 1.1.10)
# 'id' is 'chrN:start-end'
# 'promoters' contains gene name with '-promoter' suffix
# 'genes' and 'promoters' can have multiple values separated by ',' (without space) and also 'NA'
CNVs <- fread(paste0(CNVDIR, "/out/cnv-table-gt-refpool.csv"), sep="\t", data.table=FALSE, header = TRUE)
CNVs <- CNVs[, c(1:3,6,12,14)]
CNVs$promoters <- gsub("-promoter", "", CNVs$promoters)
CNVs$genes[is.na(CNVs$genes)] <- CNVs$promoters[is.na(CNVs$genes)]
CNVs <- CNVs[unlist(lapply(CNVs$genes, function(x) gene %in% unlist(strsplit(x, split=",")))) & CNVs$sample == sample,]
if (nrow(CNVs) > 0) {
    CNVs$label <- paste0(CNVs$type, " BF=", CNVs$BF)
    CNVs[,c("start", "end")] <- t(sapply(CNVs$id, FUN = function(x) as.integer(unlist(strsplit(unlist(strsplit(x, ":"))[2], "-")))))
    print(CNVs)
    ## Getting a label for a multi-gene CNV:
    for (i in 1:nrow(CNVs)) {
        if ((CNVs[i, "start"] < window.start) || (CNVs[i, "end"] > window.end)) {
            CNVs[i, "label"] <- paste0("multi-gene ", CNVs[i, "label"])
            CNVs[i, "start"] <- max(CNVs[i, "start"], window.start)
            CNVs[i, "end"] <- min(CNVs[i, "end"], window.end)
        }
    }
}

## The gene's panel actual capture
captures <- fread(paste0(PreProcessDIR, "/", PROJECT, "-capture.bed"), sep="\t", data.table=FALSE, header = FALSE, col.names = df.header[1:3])
captures$start <- captures$start + 1

## ------------------------------------------------------------------------
## Relative coverage

## Gene types
source("gene_type.R")

## Pedigree and reference pools information
# The header of the file is:
# ped  person  gender  pool  sexspecpool
# 'ped' is the family ID: sample ID for an unrelated sample or ID of one of the related samples
# 'person' is the sample ID
# 'gender' is the inferred sample's sex
# 'pool' is the list of unrelated samples to choose the reference set for CNV calling from; ','-separetad (no space)
# 'sexspecpool' is the list of same-sex unrelated samples to choose the reference set for CNV calling from; ','-separetad (no space)
sample.info <- read.delim(paste0(CNVDIR, "/ped-refpools.tsv"), as.is=TRUE)
sample.gender <- sample.info[sample.info$person==sample, "gender"]

## Difference between autosomal+PAR and sex genes:
if (gene %in% auto.genes) {
    pool.name <- "pool"
    table.pool.name <- "Autosomal.CNV.reference.set"
    mean.name <- "auto.genes"
} else {
    print("Gene on X")
    pool.name <- "sexspecpool"
    table.pool.name <- "X.Chr.CNV.reference.set"
    mean.name <- "sex.genes"
}

## Mean exonic coverage for each sample (created by mean_coverage.R)
mean.cov.filename <- paste0(CNVDIR, "/", BATCH, "-mean-cov.rda")
cat("Loading mean coverage for this batch... ", file=stdout())
load(mean.cov.filename)
cat("Done\n", file=stdout())
sample.exonic.mean <- sample.mean.cov[sample,mean.name]
sample.exonic.mean.percentile <- as.integer(cut(sample.exonic.mean, quantile(sample.mean.cov[, mean.name], probs=0:100/100), include.lowest=TRUE))

## CNVs summary info
# The fields we need from this table are:
# Sample.ID  Autosomal.CNV.reference.set  X.Chr.CNV.reference.set
# 'Autosomal.CNV.reference.set' is 10 samples selected as reference set for CNV calling in the autosomal genes
# 'X.Chr.CNV.reference.set' is 10 samples selected as reference set for CNV calling in the X-chromosome genes
CNVs.info <- read.delim(paste0(CNVDIR, "/out/CNVsummary.tsv"))

## Gene coverage
# 
cov <- fread(paste0(covDIR, "/merge/gene/", gene, ".cov"), sep=" ", fill=TRUE, data.table=FALSE, header=TRUE)
# The input file for each gene contains a table with coverage of all gene's bases (including intronic) with the flanks for each sample.
# Structure of this file (space-separated):
# loc                Sample1 Sample2 Sample3 Sample4 ...
# start_gene-flank   112     94      171     140     ...
# start_gene-flank+1 114     96      175     145     ...
# ...
# end_gene+flank     119     98      185     154     ...
print(dim(cov))

## Reference samples
if (gene %in% auto.genes) {
    ref.samples.all <- colnames(cov)[-match(c("loc", sample), colnames(cov))]
} else {
    ref.samples.all <- sample.info[sample.info$gender==sample.gender & sample.info$person != sample, "person"]
}
ref.samples.pool <- unlist(strsplit(sample.info[sample.info$person==sample, pool.name], ","))
ref.samples.set <- unlist(strsplit(unlist(strsplit(as.character(CNVs.info[match(sample, CNVs.info[, "Sample.ID"]), table.pool.name]), " "))[1], ","))
print(length(ref.samples.all))
print(length(ref.samples.pool))
print(length(ref.samples.set))

## This gene in CNV-calling regions
indexes <- c()
for (i in 1:length(CNVregions$start)) {
    indexes <- c(indexes, match(CNVregions$start[i], cov$loc):match(CNVregions$end[i], cov$loc))
}
cov$reg.type <- "nonCNV"
cov$reg.type[indexes] <- "CNV"
cov <- cov[,c("loc","reg.type", sample, ref.samples.all)]
print(dim(cov))

## Coordinates of (intronic) regions where we don't call CNVs:
empty.loc <- data.frame(start=c(cov[1,1], CNVregions$end + 1), end=c(CNVregions$start - 1, cov[nrow(cov),1]))
empty.loc <- empty.loc[empty.loc$end - empty.loc$start >= 0,]
## testing if it is correct:
e.indexes <- c()
for (i in 1:length(empty.loc$start)) {
    e.indexes <- c(e.indexes, match(empty.loc$start[i], cov$loc):match(empty.loc$end[i], cov$loc))
}
if (! identical(sort(c(indexes, e.indexes)), 1:nrow(cov))) {
    stop("Coordinates of CNV and nonCNV regions do not match the whole")
}

## Relative coverage
relative.cov <- cov[cov$reg.type=="CNV",]
for (s in c(sample, ref.samples.all)) {
    relative.cov[, s] <- relative.cov[, s]/sample.mean.cov[s, mean.name]
}
mean.relative.cov <- apply(relative.cov[, ref.samples.pool], 1, mean)

## Normalized coverage
norm.relative.cov <- relative.cov
for (s in c(sample, ref.samples.all)) {
    norm.relative.cov[, s] <- norm.relative.cov[, s]/mean.relative.cov
}

non.ref.samples <- setdiff(ref.samples.all, ref.samples.set)

## ------------------------------------------------------------------------
## Plot

## Load normal coverage df.gr.gene.cov.summary:
load(paste0(covDIR, "/geneplot/", gene, ".rda"))
df.gr.gene.cov.summary.soi <- cbind(df.gr.gene.cov.summary, soi=cov[, sample])

gr.gene.cov.summary <- as(df.gr.gene.cov.summary.soi, "GRanges")

sample.legend <- paste0(sample, " (mean exonic cov: ", round(sample.exonic.mean), "/percentile: ", sample.exonic.mean.percentile, ")")
dTrack.gene.summary <- DataTrack(gr.gene.cov.summary,
                                 name = "Coverage",
                                 type = c("l", "g"),
                                 groups = factor(c("Lower (5th percentile)", "Median (50th percentile)", "Upper (95th percentile)", sample.legend),
                                     c("Lower (5th percentile)", "Median (50th percentile)", "Upper (95th percentile)", sample.legend)),
                                 col = c("red", "#0099FF", "darkgreen", "black"),
                                 legend = TRUE,
                                 cex.legend = 0.8,
                                 fontsize.legend = 30,
                                 fontcolor.legend = "black",
                                 cex.axis = 1.7,
                                 lwd = c(1,1,1,2)
                                 )
track.list <- list(dTrack.gene.summary)

## Coverage problematic regions, if any:
if(nrow(pr5p) > 0) {
    atrack.pr5p <- AnnotationTrack(start = pr5p$start, end = pr5p$end,
                                   strand="*", chromosome = chr,
                                   name = "Prob. regions", fill="red", stacking = "dense")
    track.list <- c(track.list, atrack.pr5p)
}

## Target and capture
atrack.target <- AnnotationTrack(start = targets$start, end = targets$end,
                                 strand="*", chromosome = chr,
                                 name = "Target",  stacking = "dense")
atrack.capture <- AnnotationTrack(start = captures$start, end = captures$end,
                                  strand="*", chromosome = chr,
                                  name = "Capture",  stacking = "dense", fill = "lightgreen")
track.list <- c(track.list, atrack.target, atrack.capture)

## Normalized coverage in CNV calling regions
df.gene.cov <- cbind(data.frame(chr=chr, start=norm.relative.cov$loc, end=norm.relative.cov$loc, strand="*"),
    norm.relative.cov[,c(non.ref.samples, ref.samples.set, sample)])
## Plus intronic regions with NA:
intronic <- cbind(data.frame(chr=chr, start=empty.loc$start, end=empty.loc$end, strand="*"),
                  data.frame(matrix(NA, nrow(empty.loc), length(c(non.ref.samples, ref.samples.set, sample)))))
colnames(intronic)[-(1:4)] <- c(non.ref.samples, ref.samples.set, sample)
df.gene.cov <- rbind(df.gene.cov, intronic)
df.gene.cov <- df.gene.cov[order(df.gene.cov$start),]
gr.gene.cov <- as(df.gene.cov, "GRanges")

dTrack.gene <- DataTrack(gr.gene.cov,
                         name = "Normalized relative coverage",
                         type = c("l", "g"),
                         groups = factor(c(rep("non-reference samples", length(non.ref.samples)), rep("reference samples", length(ref.samples.set)), sample),
                             levels = c("non-reference samples", "reference samples", sample)),
                         col=c("gray60", "#3399CC", "black"),
                         ## groups = factor(c(non.ref.samples, ref.samples.set, sample),
                         ##     levels=c(non.ref.samples, ref.samples.set, sample)),
                         ## col=c(rep("gray60", length(non.ref.samples)), rep("#3399CC", length(ref.samples.set)), "black"),
                         cex.legend=0.8,
                         fontsize.legend=30,
                         fontcolor.legend="black",
                         cex.axis=1.7,
                         ylim=c(0,2)
                         )
track.list <- c(track.list, dTrack.gene)

## Problematic CNV calling regions, if any:
if(nrow(CNVprregions) > 0) {
    atrack.CNVprregions <- AnnotationTrack(start = CNVprregions$start, end = CNVprregions$end,
                                       strand="*", chromosome = chr,
                                       name = "CNV prob. regions", fill="red", stacking = "dense")
    track.list <- c(track.list, atrack.CNVprregions)
}

## CNV calling regions
atrack.CNVregions <- AnnotationTrack(start = CNVregions$start, end = CNVregions$end,
                                     strand="*", chromosome = chr,
                                     name = "CNV calling regions", fill="darkgreen", stacking = "dense")
track.list = c(track.list, atrack.CNVregions)

## CNVs in this sample
atrack.CNVs <- AnnotationTrack(name = "sample CNVs", start = CNVs$start, end = CNVs$end, strand="*", chromosome = chr,
                               group = CNVs$label, stacking = "dense", id = CNVs$label)
displayPars(atrack.CNVs) <- list(showFeatureId=TRUE, fontcolor.item="black", fontsize.item=20)
track.list <- c(track.list, atrack.CNVs)

## Transcript
ensembl.grch37.mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",
                               host="grch37.ensembl.org",
                               path="/biomart/martservice",
                               dataset="hsapiens_gene_ensembl")
fm <- Gviz:::.getBMFeatureMap()
fm["symbol"] <- "external_gene_name"
biomTrack <- BiomartGeneRegionTrack(genome=ref, chromosome=chr, start=window.start, end=window.end, name=gene,
                                    featureMap=fm, filters=list(ensembl_transcript_id = transcript.id),
                                    biomart=ensembl.grch37.mart, transcriptAnnotation = "transcript",
                                    just.group = "above", fontsize.group=30, fontcolor.group="black", fontsize.title=12)
track.list <- c(track.list, biomTrack)

itrack <- IdeogramTrack(genome=ref, chromosome=chr, fontsize=30, fontcolor="black")
gtrack <- GenomeAxisTrack(fontsize=20, fontcolor="black")

track.list <- c(track.list, gtrack, itrack)

h <- 20
# plot width is propotional to gene length, set a factor, p.fac=1500 based on length of ANXA5 (~=30000) / width (20)
p.fac <- 1500
w <- (window.end - window.start)/p.fac
if (w < 26) {
    w <- 26
} else if (w > 200) {
  w <- 200
}

png(file = paste0(CNVDIR, "/out/CNVpics/", sample, "-", gene, ".png"), width = w, height = h, units="in", res=100)
plotTracks(track.list, frame=FALSE, from = window.start, to = window.end,
           col.axis="black", col.title="black", background.title="lightgray",
           cex.title=1.5)
dev.off()
