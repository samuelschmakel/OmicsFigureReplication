library(ggplot2)
library(DESeq2)
library(GEOquery)

# Load file (separator is space by default)
counts <- as.matrix(read.delim2("data/GSE282850_tovar-nishino_rnaseq_raw_counts.txt"))

# Get the metadat from GEO
gse <- getGEO("GSE282850")  # Replace with your GSE ID
meta <- pData(phenoData(gse[[1]]))
