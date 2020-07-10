## ---- include=FALSE-----------------------------------------------------------
library(knitr)
opts_chunk$set(out.width="100%", cache=TRUE)

## ----setup, echo=FALSE--------------------------------------------------------
suppressPackageStartupMessages({ 
    library(CNVRanger)
    library(AnnotationHub)
    library(regioneR)
    library(BSgenome.Btaurus.UCSC.bosTau6.masked)
    library(SummarizedExperiment)
    library(curatedTCGAData)
    library(TCGAutils)
    library(RaggedExperiment)
})

## -----------------------------------------------------------------------------
showClass("RaggedExperiment")

## -----------------------------------------------------------------------------
sample1 <- GRanges(
    c(A = "chr1:1-10:-", B = "chr1:8-14:+", C = "chr1:15-18:+"),
    score = 3:5, type=c("germline", "somatic", "germline"))
sample2 <- GRanges(
    c(D = "chr1:1-10:-", E = "chr1:11-18:+"),
    score = 11:12, type=c("germline", "somatic"))

## -----------------------------------------------------------------------------
colDat <- DataFrame(id=1:2, status = factor(c("control", "case")))

## -----------------------------------------------------------------------------
(ragexp <- RaggedExperiment(
    sample1 = sample1,
    sample2 = sample2,
    colData = colDat))

## -----------------------------------------------------------------------------
grl <- GRangesList(sample1=sample1, sample2=sample2)
ragexp2 <- RaggedExperiment(grl, colData = colDat)
identical(ragexp, ragexp2)

## -----------------------------------------------------------------------------
rowRanges(ragexp)

## ---- echo=FALSE, fig.cap="RaggedExperiment object schematic. Rows and columns represent genomic ranges and samples, respectively. Assay operations can be performed with (from left to right) compactAssay, qreduceAssay, and sparseAssay.", out.width="\\maxwidth"----
knitr::include_graphics("RaggedExperiment.png")

## -----------------------------------------------------------------------------
sparseAssay(ragexp)

## -----------------------------------------------------------------------------
unlist(grl)

## -----------------------------------------------------------------------------
assay(ragexp, "score")
assay(ragexp, "type")

## -----------------------------------------------------------------------------
compactAssay(ragexp)
compactAssay(ragexp, "type")

## -----------------------------------------------------------------------------
disjoinAssay(ragexp, simplifyDisjoin = mean)

## -----------------------------------------------------------------------------
weightedmean <- function(scores, ranges, qranges)
{
    isects <- pintersect(ranges, qranges)
    sum(scores * width(isects)) / sum(width(isects))
}

## -----------------------------------------------------------------------------
grl

## -----------------------------------------------------------------------------
(query <- GRanges(c("chr1:1-14:-", "chr1:15-18:+")))

## -----------------------------------------------------------------------------
qreduceAssay(ragexp, query, simplifyReduce = weightedmean)

## ---- results='hide'----------------------------------------------------------
sparseSummarizedExperiment(ragexp)
compactSummarizedExperiment(ragexp)
disjoinSummarizedExperiment(ragexp, simplify = mean)
qreduceSummarizedExperiment(ragexp, query = query, simplify = weightedmean)

## ---- echo=FALSE--------------------------------------------------------------
knitr::include_graphics("CNVRanger.png")

## ----input, echo=FALSE, fig.wide=TRUE-----------------------------------------
vign.dir <- system.file("vignettes", package = "CNVRanger")
knitr::include_graphics(file.path(vign.dir, "Input.png"))

## ----lib----------------------------------------------------------------------
library(CNVRanger)

## ----readCalls----------------------------------------------------------------
data.dir <- system.file("extdata", package="CNVRanger")
call.file <- file.path(data.dir, "Silva16_PONE_CNV_calls.csv")
calls <- read.csv(call.file, as.is=TRUE)
nrow(calls)
head(calls)

## ----nrSamples----------------------------------------------------------------
length(unique(calls[,"NE_id"]))

## ----cnvCalls-----------------------------------------------------------------
grl <- makeGRangesListFromDataFrame(calls, 
    split.field = "NE_id", keep.extra.columns = TRUE)
grl

## ----sortCalls----------------------------------------------------------------
grl <- sort(grl)
grl

## ----RaggedExperiment---------------------------------------------------------
ra <- RaggedExperiment(grl)
ra

## ----RaggedExperiment-assay---------------------------------------------------
assay(ra[1:5,1:5])

## ----RaggedExperiment-colData-------------------------------------------------
weight <- rnorm(ncol(ra), mean=1100, sd=100)
fcr <- rnorm(ncol(ra), mean=6, sd=1.5)
colData(ra)$weight <- round(weight, digits=2)
colData(ra)$fcr <- round(fcr, digits=2)
colData(ra)

## ----summarization, echo=FALSE, out.width = "40%", out.extra='style="float:right; padding:10px"'----
vign.dir <- system.file("vignettes", package = "CNVRanger")
knitr::include_graphics(file.path(vign.dir, "Summarization.png"))

## ----cnvrs--------------------------------------------------------------------
cnvrs <- populationRanges(grl, density = 0.1)
cnvrs

## ----gistic-------------------------------------------------------------------
cnvrs <- populationRanges(grl, density = 0.1, est.recur = TRUE)
cnvrs

## ----recurr-------------------------------------------------------------------
subset(cnvrs, pvalue < 0.05)

## ----recurrViz----------------------------------------------------------------
plotRecurrentRegions(cnvrs, genome = "bosTau6", chr = "chr1")

## ----overlap, echo=FALSE, out.width = "40%", out.extra='style="float:right; padding:10px"'----
vign.dir <- system.file("vignettes", package = "CNVRanger")
knitr::include_graphics(file.path(vign.dir, "Overlap.png"))

## ----getBtGenes---------------------------------------------------------------
library(AnnotationHub)
ah <- AnnotationHub::AnnotationHub()
ahDb <- AnnotationHub::query(ah, pattern = c("Bos taurus", "EnsDb"))
ahDb

## ----getBtGenes2--------------------------------------------------------------
ahEdb <- ahDb[["AH60948"]]
bt.genes <- ensembldb::genes(ahEdb)
GenomeInfoDb::seqlevelsStyle(bt.genes) <- "UCSC"
bt.genes

## ----formatBtGenes------------------------------------------------------------
sel.genes <- subset(bt.genes, seqnames %in% paste0("chr", 1:2))
sel.genes <- subset(sel.genes, gene_biotype == "protein_coding")
sel.cnvrs <- subset(cnvrs, seqnames %in% paste0("chr", 1:2))

## ----findOlaps----------------------------------------------------------------
olaps <- GenomicRanges::findOverlaps(sel.genes, sel.cnvrs, ignore.strand=TRUE)
qh <- S4Vectors::queryHits(olaps)
sh <- S4Vectors::subjectHits(olaps)
cgenes <- sel.genes[qh]
cgenes$type <- sel.cnvrs$type[sh]
subset(cgenes, select = "type")

## ----vizOlaps-----------------------------------------------------------------
cnvOncoPrint(grl, cgenes)

## ----ovlpTest-----------------------------------------------------------------
library(regioneR)
library(BSgenome.Btaurus.UCSC.bosTau6.masked)
res <- suppressWarnings(overlapPermTest(A=sel.cnvrs, B=sel.genes, ntimes=100, 
    genome="bosTau6", mask=NA, per.chromosome=TRUE, count.once=TRUE))
res

## ----permDist-----------------------------------------------------------------
summary(res[[1]]$permuted)

## ----vizPermTest--------------------------------------------------------------
plot(res)

## ----expression, echo=FALSE, out.width = "40%", out.extra='style="float:right; padding:10px"'----
vign.dir <- system.file("vignettes", package = "CNVRanger")
knitr::include_graphics(file.path(vign.dir, "Expression.png"))

## ----rseqdata-----------------------------------------------------------------
rseq.file <- file.path(data.dir, "counts_cleaned.txt")
rcounts <- read.delim(rseq.file, row.names=1, stringsAsFactors=FALSE)
rcounts <- as.matrix(rcounts)
dim(rcounts)
rcounts[1:5, 1:5]

## ----rse----------------------------------------------------------------------
library(SummarizedExperiment)
rse <- SummarizedExperiment(assays=list(rcounts=rcounts), 
                rowRanges=granges(sel.genes)[rownames(rcounts)])
rse

## ----cnvEQTL------------------------------------------------------------------
res <- cnvEQTL(cnvrs, grl, rse, window="1Mbp", verbose=TRUE)
res

## ----tcgaSetup, message = FALSE-----------------------------------------------
library(curatedTCGAData)
gbm <- curatedTCGAData::curatedTCGAData("GBM",
        assays=c("GISTIC_Peaks", "CNVSNP", "RNASeq2GeneNorm"),
        dry.run=FALSE)
gbm

## ----tcgaGeneAnno, message = FALSE--------------------------------------------
library(TCGAutils)
gbm <- TCGAutils::symbolsToRanges(gbm, unmapped=FALSE)
for(i in 1:3)
{
    rr <- rowRanges(gbm[[i]])
    GenomeInfoDb::genome(rr) <- "NCBI37"
    GenomeInfoDb::seqlevelsStyle(rr) <- "UCSC"
    ind <- as.character(seqnames(rr)) %in% c("chr1", "chr2")
    rowRanges(gbm[[i]]) <- rr
    gbm[[i]] <- gbm[[i]][ind,]
}
gbm

## ----gbmIntersect-------------------------------------------------------------
gbm <- MultiAssayExperiment::intersectColumns(gbm)
TCGAutils::sampleTables(gbm)
data(sampleTypes, package="TCGAutils")
sampleTypes
gbm <- TCGAutils::splitAssays(gbm, sampleCodes="01")
gbm

## ----segmean------------------------------------------------------------------
ind <- grep("CNVSNP", names(gbm))
head( mcols(gbm[[ind]]) )
summary( mcols(gbm[[ind]])$Segment_Mean )

## ----lr2int, eval=FALSE-------------------------------------------------------
#  round( (2^lr) * 2)

## ----transformToStates--------------------------------------------------------
smean <- mcols(gbm[[ind]])$Segment_Mean
state <- round(2^smean * 2)
state[state > 4] <- 4
mcols(gbm[[ind]])$state <- state
gbm[[ind]] <- gbm[[ind]][state != 2,]
mcols(gbm[[ind]]) <- mcols(gbm[[ind]])[,3:1]
table(mcols(gbm[[ind]])$state)

## ----gbmEQTL------------------------------------------------------------------
res <- cnvEQTL(cnvrs = "01_GBM_GISTIC_Peaks-20160128", 
               calls = "01_GBM_CNVSNP-20160128", 
               rcounts = "01_GBM_RNASeq2GeneNorm-20160128_ranged", 
               data = gbm, window = "1Mbp", verbose = TRUE)
res

## ----plotEQTL-----------------------------------------------------------------
res[2]
(r <- GRanges(names(res)[2]))
plotEQTL(cnvr=r, genes=res[[2]], genome="hg19", cn="CN1")

## ----echo=FALSE, fig.height=2, fig.width=2------------------------------------
knitr::include_graphics("overview_1.png")

## ----eval=FALSE---------------------------------------------------------------
#  java -jar mutect.jar \
#      --analysis_type MuTect \
#      -R hg38.fasta \
#      --dbsnp $DBSNP_VCF \
#      --cosmic $COSMIC_VCF \
#      -I:tumor $BAM_TUMOR  \
#      -o $OUT/${SAMPLEID}_mutect_stats.txt \
#      -vcf $OUT/${SAMPLEID}_mutect.vcf

## ----eval=FALSE---------------------------------------------------------------
#  Rscript $PURECN/PureCN.R \
#      --out $PCN_OUT/PureCN/tumor_only/$SAMPLEID \
#      --tumor $PCN_OUT/tumor_cov/${SAMPLEID}_coverage_loess.txt \
#      --SAMPLEID ${SAMPLEID} \
#      --vcf $MUTECT_OUT/stat_tumor_only/${SAMPLEID}_mutect.vcf \
#      --statsfile $MUTECT_OUT/stat_tumor_only/${SAMPLEID}_mutect_stats.txt \
#      --normaldb $PCN_OUT/normalDB/normalDB_hg38.rds \
#      --normal_panel $PCN_OUT/normalDB/mapping_bias_hg38.rds \
#      --intervals $INPUT/bed/baits_hg38_intervals.txt \
#      --intervalweightfile $PCN_OUT/normalDB/interval_weights_hg38.txt \
#      --snpblacklist hg38_simpleRepeats.bed \
#      --genome hg38 \
#      --force --postoptimize --seed 123

## ----echo=FALSE, out.width = '180%'-------------------------------------------
knitr::include_graphics("overview_2.png")

## ----echo=FALSE, out.width = '180%'-------------------------------------------
knitr::include_graphics("overview_3.png")

## ----example_output, echo=FALSE, warning=FALSE, message=FALSE-----------------
library(PureCN)
file.rds = "Sample1_PureCN.rds"
ret = readRDS(file.rds)

## ----fig.width=4.5, fig.height=4.5--------------------------------------------
plotAbs(ret, type="overview")

## ----fig.width=4, fig.height=4------------------------------------------------
plotAbs(ret, 1, type="hist")

## -----------------------------------------------------------------------------
createCurationFile(file.rds)
read.csv("Sample1_PureCN.csv")

## -----------------------------------------------------------------------------
gene.calls = callAlterations(ret)
head(gene.calls)

## -----------------------------------------------------------------------------
loh = callLOH(ret)
head(loh)

## ----analysis_example, echo=FALSE---------------------------------------------
knitr::include_graphics("analysis_example.png")

