library(ggplot2)
library(DESeq2)
library(GEOquery)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)

rm(list = ls(all.names = TRUE))
gc()

source("scripts/helper_functions.R")

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

# Generate results
alpha <- 0.05
lfcThreshold <- 1

res_diff_vs_undiff <- results(dds, contrast = c("condition", "diff", "undiff"), alpha = alpha, lfcThreshold = lfcThreshold)
res_AICAR_vs_undiff <- results(dds, contrast = c("condition", "AICAR", "undiff"), alpha = alpha, lfcThreshold = lfcThreshold)
res_palmitate_vs_undiff <- results(dds, contrast = c("condition", "palmitate", "undiff"), alpha = alpha, lfcThreshold = lfcThreshold)
res_diff_vs_AICAR <- results(dds, contrast = c("condition", "diff", "AICAR"), alpha = alpha, lfcThreshold = lfcThreshold)
res_diff_vs_palmitate <- results(dds, contrast = c("condition", "diff", "palmitate"), alpha = alpha, lfcThreshold = lfcThreshold)
res_AICAR_vs_palmitate <- results(dds, contrast = c("condition", "AICAR", "palmitate"), alpha = alpha, lfcThreshold = lfcThreshold)

#### PCA Plot ####
normalized_counts <- PCA_plotting(dds)

#### DEG list output ####
diff_vs_undiff_DEG_counts <- get_counts(res_diff_vs_undiff, normalized_counts, alpha, lfcThreshold)
AICAR_vs_undiff_DEG_counts <- get_counts(res_AICAR_vs_undiff, normalized_counts, alpha, lfcThreshold)
palmitate_vs_undiff_DEG_counts <- get_counts(res_palmitate_vs_undiff, normalized_counts, alpha, lfcThreshold)
diff_vs_AICAR_DEG_counts <- get_counts(res_diff_vs_AICAR, normalized_counts, alpha, lfcThreshold)
diff_vs_palmitate_DEG_counts <- get_counts(res_AICAR_vs_undiff, normalized_counts, alpha, lfcThreshold)
AICAR_vs_palmitate_DEG_counts <- get_counts(res_AICAR_vs_palmitate, normalized_counts, alpha, lfcThreshold)

# Only keep samples that are involved in this comparison
diff_vs_undiff_DEG_counts <- diff_vs_undiff_DEG_counts[,1:8]
AICAR_vs_undiff_DEG_counts <- AICAR_vs_undiff_DEG_counts[,c(1:4,9:11)]
palmitate_vs_undiff_DEG_counts <- palmitate_vs_undiff_DEG_counts[,c(1:4,12:15)]
diff_vs_AICAR_DEG_counts <- diff_vs_AICAR_DEG_counts[,5:11]
diff_vs_palmitate_DEG_counts <- diff_vs_palmitate_DEG_counts[,c(5:8,12:15)]
AICAR_vs_palmitate_DEG_counts <- AICAR_vs_palmitate_DEG_counts[,9:15]

# Get logfoldchange information from the res object for each comparison:
diff_vs_undiff_DEGs <- as.data.frame(res_diff_vs_undiff)

# Now, set up the counts matrix for the heatmap, and create the annotation row (saved in the returned list):
diff_vs_undiff_heatmap_list <- reorder_counts(diff_vs_undiff_DEG_counts, diff_vs_undiff_DEGs)

# Create heatmap
pheatmap(
  diff_vs_undiff_heatmap_list$cts, 
  cluster_rows = TRUE, 
  cluster_cols = TRUE, 
  annotation_row = diff_vs_undiff_heatmap_list$annotation,
  scale = "row",  # Normalize rows to Z-scores
  color = colorRampPalette(c("blue", "white", "red"))(50),
  main = "Heatmap of DEGs"
)

#### GO Enrichment Analysis ####
#GO_enrichment(diff_vs_undiff_DEG_counts, "diff_vs_undiff")
