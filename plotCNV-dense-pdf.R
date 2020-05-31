library(data.table)
library(ggplot2)
library(Biobase)
library(biomaRt)
library(Gviz)
library(GenomicRanges)

config <- "/path/to/conf.R"
gene <- "G1"
sample <- "S1"

source(config)

## ------------------------------------------------------------------------
## Tracks
ref <- "hg19"
df.header <- c("chr", "start", "end", "symbol")

## This gene target and boundaries
targets <- read.delim(paste0(PreProcessDIR, "/", PROJECT, "_target_mapped_genes.filtered.bed"), header=F, col.names = c("chrom", "start", "end", "annotation"))
genetargets <- targets[sapply(targets$annotation, function(x) gene %in% unlist(strsplit(as.character(x), split=";"))),]
genetargets$start <- genetargets$start + 1
chr <- paste0("chr", genetargets[1,1])
print(chr)

design <- read.delim(paste0(PreProcessDIR, "/", PROJECT, "-genes-design.bed"), header=F, col.names = c("chrom", "start", "end", "gene"))
window.start <- design[design$gene == gene, "start"] + 1
window.end <- design[design$gene == gene, "end"]

## Problematic regions:
pr5p <- read.delim(paste0(covDIR, "/", PROJECT, "-problematic.5p.interval"), header=F, col.names = df.header)
genepr5p <- pr5p[pr5p$symbol==gene,]

## Ensembl transcript ID
transcripts <- read.delim(paste0(PreProcessDIR, "/", PROJECT, "-transcript-gene.info"), header=F, col.names = c("transcriptid", "symbol"))
transcript.id <- transcripts[transcripts$symbol == gene,]$transcriptid

## CNV calling regions
CNVregions <- fread(paste0(CNVDIR, "/cnv-regions.bed"), sep="\t", data.table=FALSE, header=FALSE, col.names=c("chrom", "start", "end", "annotation"))
CNVregions <- CNVregions[grep(paste0("^", gene, "-"), CNVregions$annotation),]
CNVregions$start <- CNVregions$start + 1
print(CNVregions)

## Exons
Exons <- fread(paste0(CNVDIR, "/cnv-exons.bed"), sep="\t", data.table=FALSE, header=FALSE, col.names=c("chrom", "start", "end", "annotation"))
Exons <- Exons[grep(paste0("^", gene, "-"), Exons$annotation),]
Exons$start <- Exons$start + 1
Exons$num <- as.numeric(sub(paste0("^", gene, "-"), "", Exons$annotation))
print(Exons)
if (Exons$num[1] == 1) {
    strand = "+"
} else {
    strand = "-"
}
print(strand)

## CNV problematic regions:
CNVpr <- read.delim(paste0(CNVDIR, "/out/poorly_covered_regions.tsv"), header = TRUE)
CNVprregions <- merge(CNVpr, data.frame(space=CNVregions$chr, start=CNVregions$start, end=CNVregions$end))
print(CNVprregions)

## Sample's CNV in this gene:
CNVs <- read.delim(paste0(CNVDIR, "/out/cnv-table-gt-refpool.csv"), header = TRUE)
CNVs <- CNVs[, c(1:3,6,12,14)]
colnames(CNVs)[6] <- "promoters"
CNVs$promoters <- gsub("-promoter", "", as.character(CNVs$promoters))
colnames(CNVs)[5] <- "genes"
CNVs$genes <- as.character(CNVs$genes)
CNVs$genes[is.na(CNVs$genes)] <- CNVs$promoters[is.na(CNVs$genes)]
CNVs <- CNVs[unlist(lapply(CNVs$genes, function(x) gene %in% unlist(strsplit(x, split=",")))) & CNVs$sample == sample,]
if (nrow(CNVs) > 0) {
    CNVs$label <- paste0(CNVs$type, " BF=", CNVs$BF)
    CNVs$id <- as.character(CNVs$id)
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
captures <- read.delim(paste0(PreProcessDIR, "/", PROJECT, "-capture.bed"), header = F, col.names = df.header[1:3])
captures$start <- captures$start + 1
cap.left <- min(genetargets[1,"start"], CNVregions[1, "start"])
cap.right <- max(genetargets[nrow(genetargets),"end"], CNVregions[nrow(CNVregions), "end"])
captureregions <- captures[captures$chr == as.character(genetargets[1,1]) & !((captures$end < cap.left) | (captures$start > cap.right)),]

## ------------------------------------------------------------------------
## Relative coverage

## Gene types
source("gene_type.R")

## Pedigree and reference pools information
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

## Mean exonic coverage for each sample
mean.cov.filename <- paste0(CNVDIR, "/", BATCH, "-mean-cov.rda")
cat("Loading mean coverage for this batch... ", file=stdout())
load(mean.cov.filename)
cat("Done\n", file=stdout())
sample.exonic.mean <- sample.mean.cov[sample,mean.name]
sample.exonic.mean.percentile <- as.integer(cut(sample.exonic.mean, quantile(sample.mean.cov[, mean.name], probs=0:100/100), include.lowest=TRUE))

## CNVs summary info
CNVs.info <- read.delim(paste0(CNVDIR, "/out/CNVsummary.tsv"))

## Gene coverage
cov <- fread(paste0(covDIR, "/merge/gene/", gene, ".cov"), sep=" ", fill=TRUE, data.table=FALSE, header=TRUE)
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
## Collapse large intronic regions:
gap <- 500
big.empty.loc <- empty.loc
big.empty.loc$start[2:(nrow(big.empty.loc))] <- big.empty.loc$start[2:(nrow(big.empty.loc))] + gap/2
big.empty.loc$end[1:(nrow(big.empty.loc)-1)] <- big.empty.loc$end[1:(nrow(big.empty.loc)-1)] - gap/2
big.empty.loc <- big.empty.loc[big.empty.loc$end-big.empty.loc$start+1 > 1.5*gap,]

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

itrack <- IdeogramTrack(genome=ref, chromosome=chr, fontcolor="black", showTitle=FALSE)

## Load normal coverage df.gr.gene.cov.summary:
load(paste0(covDIR, "/geneplot/", gene, ".rda"))
df.gr.gene.cov.summary.soi <- cbind(df.gr.gene.cov.summary, soi=cov[, sample])
for (i in 1:nrow(big.empty.loc)) {
    df.gr.gene.cov.summary.soi[match(big.empty.loc$start[i], df.gr.gene.cov.summary.soi$start):match(big.empty.loc$end[i], df.gr.gene.cov.summary.soi$start),c("lower", "median", "upper", "soi")] <- NA
}
print(dim(df.gr.gene.cov.summary.soi))
for (i in 1:nrow(big.empty.loc)) {
    df.gr.gene.cov.summary.soi <- df.gr.gene.cov.summary.soi[-(match(big.empty.loc$start[i]+gap/2, df.gr.gene.cov.summary.soi$start):match(big.empty.loc$end[i]-gap/2, df.gr.gene.cov.summary.soi$start)),]
}
rescale <- data.frame(loc.orig=df.gr.gene.cov.summary.soi$start)
a <- df.gr.gene.cov.summary.soi$start[1]
df.gr.gene.cov.summary.soi$start <- a:(a+nrow(df.gr.gene.cov.summary.soi) - 1)
df.gr.gene.cov.summary.soi$end <- a:(a+nrow(df.gr.gene.cov.summary.soi) - 1)
print(dim(df.gr.gene.cov.summary.soi))
print(head(df.gr.gene.cov.summary.soi))
rescale$loc.new <- df.gr.gene.cov.summary.soi$start

gr.gene.cov.summary <- as(df.gr.gene.cov.summary.soi, "GRanges")

sample.legend <- paste0(sample, " (mean cov.: ", round(sample.exonic.mean), "/percentile: ", sample.exonic.mean.percentile, ")")
dTrack.gene.summary <- DataTrack(gr.gene.cov.summary,
                                 name = "Raw coverage",
                                 type = c("l", "g"),
                                 groups = factor(c("Lower (5th percentile)", "Median (50th percentile)", "Upper (95th percentile)", sample.legend),
                                     c("Lower (5th percentile)", "Median (50th percentile)", "Upper (95th percentile)", sample.legend)),
                                 col = c("red", "#0099FF", "darkgreen", "black"),
                                 legend = TRUE,
                                 fontcolor.legend = "black",
                                 fontface.title = 1,
                                 )
track.list <- c(dTrack.gene.summary)

## Coverage problematic regions, if any:
if(nrow(genepr5p) > 0) {
    atrack.pr5p <- AnnotationTrack(start = rescale$loc.new[match(genepr5p$start, rescale$loc.orig)], end = rescale$loc.new[match(genepr5p$end, rescale$loc.orig)],
                                   strand="*", chromosome = chr, col=NULL,
                                   name = "Prob. reg.", fill="red", fontface.title = 1,)
    track.list <- c(track.list, atrack.pr5p)
}

## Exon
if(nrow(Exons) > 0) {
    atrack.exons <- AnnotationTrack(start = rescale$loc.new[match(Exons$start, rescale$loc.orig)], end = rescale$loc.new[match(Exons$end, rescale$loc.orig)],
                                    strand="*", chromosome = chr,
                                    name = gene,
                                    group = "Exons",
                                    fill="#FF9900", col=NULL,
                                    id = Exons$num,
                                    fontcolor.group="black",
                                    fontface.title = 4)
    displayPars(atrack.exons) <- list(showFeatureId=TRUE, showId=FALSE, fontcolor.item="black", showTitle=TRUE)

    ## Removed intronic bases:
    big.empty.loc.2 <- big.empty.loc[big.empty.loc$start>=Exons$start[1] & big.empty.loc$end<=Exons$end[nrow(Exons)],]
    if (nrow(big.empty.loc.2) > 0) {
        num.of.bp <- big.empty.loc.2$end - big.empty.loc.2$start + 1
        removed <- AnnotationTrack(start = rescale$loc.new[match(big.empty.loc.2$start, rescale$loc.orig)], end = rescale$loc.new[match(big.empty.loc.2$end, rescale$loc.orig)],
                                   strand="*", chromosome = chr, col=NULL,
                                   col.line = NA, id = num.of.bp,
                                   fontcolor.group = "black", showTitle = FALSE)
        displayPars(removed) <- list(showFeatureId=TRUE, showId=FALSE, fontcolor.item="black", rotation.item=90)
        ot <- OverlayTrack(trackList = list(atrack.exons, removed))
        track.list <- c(track.list, ot)
    } else {
        track.list <- c(track.list, atrack.exons)
    }
}

## Normalized coverage in CNV calling regions
non.ref.samples.max <- apply(norm.relative.cov[,non.ref.samples], 1, max)
non.ref.samples.min <- apply(norm.relative.cov[,non.ref.samples], 1, min)
ref.samples.set.max <- apply(norm.relative.cov[,ref.samples.set], 1, max)
ref.samples.set.min <- apply(norm.relative.cov[,ref.samples.set], 1, min)

df.gene.cov <- cbind(data.frame(chr=chr, start=norm.relative.cov$loc, end=norm.relative.cov$loc, strand="*"),
                     cbind(non.ref.samples.max, non.ref.samples.min, ref.samples.set.max, ref.samples.set.min, norm.relative.cov[, sample]))
colnames(df.gene.cov)[-(1:4)] <- c("non.ref.samples.max", "non.ref.samples.min", "ref.samples.set.max", "ref.samples.set.min", sample)
## Plus intronic regions with NA:
intronic <- cbind(data.frame(chr=chr, start=empty.loc$start, end=empty.loc$end, strand="*"),
                  data.frame(matrix(NA, nrow(empty.loc), 5)))
colnames(intronic)[-(1:4)] <- c("non.ref.samples.max", "non.ref.samples.min", "ref.samples.set.max", "ref.samples.set.min", sample)
df.gene.cov <- rbind(df.gene.cov, intronic)
df.gene.cov <- df.gene.cov[order(df.gene.cov$start),]
## rescale:
df.gene.cov$start <- rescale$loc.new[match(df.gene.cov$start, rescale$loc.orig)]
df.gene.cov$end <- df.gene.cov$start
gr.gene.cov <- as(df.gene.cov, "GRanges")

dTrack.gene <- DataTrack(gr.gene.cov,
                         name = "Norm. relative coverage",
                         type = c("l", "g"),
                         groups = factor(c(rep(format("Non-reference samples", width=40), 2), rep(format("Reference samples", width=40), 2), format(sample, width=40)),
                             levels = c(format("Non-reference samples", width=40), format("Reference samples", width=40), format(sample, width=40))),
                         col=c("gray60", "#3399CC", "black"),
                         legend = TRUE,
                         fontcolor.legend = "black",
                         ylim = c(0.2,1.5),
                         fontface.title = 1,
                         )
track.list <- c(track.list, dTrack.gene)

## Problematic CNV calling regions, if any:
if(nrow(CNVprregions) > 0) {
    atrack.CNVprregions <- AnnotationTrack(start = rescale$loc.new[match(CNVprregions$start, rescale$loc.orig)], end = rescale$loc.new[match(CNVprregions$end, rescale$loc.orig)],
                                           strand = "*", chromosome = chr, col = NULL,
                                           name = "CNV prob. regions", fill = "red", fontface.title = 1,)
    track.list <- c(track.list, atrack.CNVprregions)
}

## CNV calling regions
if(nrow(CNVregions) > 0) {
    atrack.CNVregions <- AnnotationTrack(start = rescale$loc.new[match(CNVregions$start, rescale$loc.orig)], end = rescale$loc.new[match(CNVregions$end, rescale$loc.orig)],
                                         strand = "*", chromosome = chr, col = "white",
                                         name = "Calling regions",
                                         fill = "springgreen3",
                                         collapse = FALSE,
                                         col.line = NA,
                                         fontcolor.group = NA,
                                         fontface.title = 1,)
    if(nrow(CNVs) > 0) {
        atrack.CNVs <- AnnotationTrack(start = rescale$loc.new[match(CNVs$start, rescale$loc.orig)], end = rescale$loc.new[match(CNVs$end, rescale$loc.orig)],
                                       strand = "*", chromosome = chr,
                                       fill = "red",
                                       id = CNVs$label,
                                       col = NULL, showTitle = FALSE)
        displayPars(atrack.CNVs) <- list(showFeatureId = TRUE, fontcolor.item = "black", fontface.item = 2, alpha.item = 1, alpha = 0.15)
        track.list <- c(track.list, OverlayTrack(trackList = list(atrack.CNVregions, atrack.CNVs)))

    } else {
        track.list <- c(track.list, atrack.CNVregions)
    }
}

## Transcript
ensembl.grch37.mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",
                               host="grch37.ensembl.org",
                               path="/biomart/martservice",
                               dataset="hsapiens_gene_ensembl")
fm <- Gviz:::.getBMFeatureMap()
fm["symbol"] <- "external_gene_name"
biomTrack <- BiomartGeneRegionTrack(genome=ref, chromosome=chr, start=window.start, end=window.end, name=gene,
                                    featureMap=fm, filters=list(ensembl_transcript_id = transcript.id),
                                    biomart=ensembl.grch37.mart, transcriptAnnotation ="transcript",
                                    min.height = 10,
                                    just.group = "above", fontcolor.group="black", showTitle=FALSE, col=NULL)
track.list.2 <- list(biomTrack)

gtrack <- GenomeAxisTrack(fontcolor="black", distFromAxis = 2, labelPos = "above")
track.list.2 <- c(track.list.2, gtrack)

track.list.2 <- c(track.list.2, itrack)

## Rescale window:
window.start.2 <- window.start
window.end.2 <- window.end
window.start <- df.gr.gene.cov.summary.soi$start[1]
window.end <- df.gr.gene.cov.summary.soi$start[nrow(df.gr.gene.cov.summary.soi)]

## A4: 8.267 Inches x 11.692 Inches
h <- 20
p.fac <- 1500
w <- (window.end - window.start)/p.fac
if (w < 42.4) {
    w <- 42.4
}
h <- h/2
w <- w/2
fontscale <- 3.3/2

pdf(file = paste0(CNVDIR, "/out/CNVpics/new-dense-", sample, "-", gene, ".pdf"), width = w, height = h)

grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 1, heights = unit(c(16*h/20,4*h/20), "in"))))

pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 1))
plotTracks(track.list,
           from = window.start, to = window.end,
           col.axis="black", col.title="black",
           lwd = 3,
           cex = fontscale,
           cex.legend = 0.85*fontscale,
           cex.title = 0.9*fontscale,
           cex.group = fontscale,
           cex.axis = 0.9*fontscale,
           sizes = c(0.4,0.12,0.36,0.12),
           title.width = 0.8,
           add = TRUE)
popViewport(1)

pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 2))
plotTracks(track.list.2,
           from = window.start.2, to = window.end.2,
           lwd = 3,
           cex = fontscale,
           cex.legend = fontscale,
           cex.title = fontscale,
           cex.group = 0.8*fontscale,
           fontface.group = 1,
           cex.axis = fontscale,
           sizes = c(0.4,0.4,0.2),
           col.axis = "black",
           add = TRUE)
popViewport(1)

dev.off()
