library(ggplot2)
library(DESeq2)
library(GEOquery)

rm(list = ls(all.names = TRUE))
gc()

# Load file (separator is space by default)
raw_data <- as.matrix(read.delim2("data/GSE282850_tovar-nishino_rnaseq_raw_counts.txt", header = TRUE, stringsAsFactors = FALSE))

cts <- apply(raw_data[, 2:ncol(raw_data)], 2, as.numeric)
rownames(cts) <- raw_data[,1]

# Get the metadata from GEO
gse <- getGEO("GSE282850")  # Replace with your GSE ID
colData <- pData(phenoData(gse[[1]]))

# Rename metadata to match counts matrix
rownames(colData) <- gsub('_mRNA', '', colData$description.1)

# Add condition column to column data
colData$condition <- factor(c("undiff","undiff","undiff","undiff","diff","diff","diff","diff","AICAR","AICAR","AICAR","palmitate","palmitate","palmitate","palmitate"), levels = c("undiff", "diff", "AICAR", "palmitate"))
colData <- colData[, c("condition"), drop = FALSE]

# Validate alignment of counts and column data
if (!all(colnames(cts) %in% rownames(colData))) stop("Sample names in countData and colData do not match.")

dds <- DESeqDataSetFromMatrix(countData = cts, colData = colData, design = ~ condition)
dds

# Pre-filtering
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

# Defining factor levels
dds$condition <- relevel(dds$condition, ref = "undiff")

# DEG analysis
dds <- DESeq(dds)
res <- results(dds, name="condition_diff_vs_undiff", alpha=0.05, lfcThreshold=.6)
res

# Calculate p-values and adjusted p-values
resOrdered <- res[order(res$pvalue),]

summary(res)

# MA-plot
plotMA(res, ylim=c(-2,2))

#### PCA Plot ####

# Perform rlog transformation
rld <- rlog(dds, blind = TRUE)


plotPCA(rld, intgroup = "condition")


# Extract transformed counts
normalized_counts <- assay(rld)

# Alternatively, use vst:
# vst <- vst(dds, blind = TRUE)
# plotPCA(vst, intgroup = "condition")
# pca_data <- assay(vst)

pca <- prcomp(t(normalized_counts))

percentVar <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 1)

# Create a data frame for ggplot
pca_df <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  condition = colData(dds)$condition
)

# Figure 2a
ggplot(pca_df, aes(x = PC1, y = PC2, color = condition, shape = condition)) +
  geom_point(size = 3, color = "black", aes(fill = condition)) +
  scale_shape_manual(values = c(24, 21, 21, 21)) +
  scale_fill_manual(values = c("red", "yellow", "green", "blue"))  +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA of Gene Expression by Condition") +
  theme_bw()

### DEG list output ###

DEGs <- as.data.frame(res[res$padj < 0.05 & abs(res$log2FoldChange) > 0.6,])
DEGsup <- DEGs[DEGs$log2FoldChange > 0,]
DEGsdown <- DEGs[DEGs$log2FoldChange < 0,]
