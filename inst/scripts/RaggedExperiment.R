library(RaggedExperiment)

## GRanges constructor function
sample1 <- GRanges(
    c(A = "chr1:1-10:-", B = "chr1:8-14:+", C = "chr1:15-18:+"),
    score = c(3, 4, 5),
    type = c("germline", "somatic", "germline")
)
sample1

## mock read-in of data
as.data.frame(sample1)
makeGRangesFromDataFrame(as.data.frame(sample1), keep.extra.columns = TRUE)

## coercion from character vector
sample2 <- as(c(A = "chr1:1-10:-", B = "chr1:11-18:+"), "GRanges")
sample2

## adding mcols data from S4Vectors DataFrame
mcols(sample2) <- DataFrame(score = c(3, 4), type = c("germline", "somatic"))

sample2

colDat <- DataFrame(id=1:2, status = factor(c("control", "case")))


ragexp <- RaggedExperiment(
    sample1 = sample1,
    sample2 = sample2,
    colData = colDat
)
ragexp

