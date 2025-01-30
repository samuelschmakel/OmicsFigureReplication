if (!requireNamespace("librarian", quietly = TRUE)) install.packages("librarian")
librarian::shelf(ggplot2, DESeq2, GEOquery, pheatmap, clusterProfiler, org.Hs.eg.db, 
                 tidyverse, ComplexHeatmap, EnhancedVolcano, RColorBrewer, gprofiler2)

rm(list = ls(all.names = TRUE))
gc()

source("scripts/helper_functions.R")

#### Setup ####

# Define directories
in_dir <- file.path("data")
work_dir <- file.path("results")
fig_dir <- file.path(work_dir, "figs")
out_dir <- file.path(work_dir, "output")

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Define file paths
raw_count_file <- file.path(in_dir, "GSE282850_tovar-nishino_rnaseq_raw_counts.txt")
metadata_file <- file.path(out_dir, "GSE282850_metadata.rds")

# Load raw counts
raw_counts <- read.delim2(file.path(in_dir, "GSE282850_tovar-nishino_rnaseq_raw_counts.txt"), row.names = 1, check.names = FALSE)

# Get metadata if not already available
if (!file.exists(metadata_file)) {
  gse <- getGEO("GSE282850", GSEMatrix = TRUE)
  colData <- pData(phenoData(gse[[1]]))
  saveRDS(colData, metadata_file)
} else {
  colData <- readRDS(metadata_file)
}

#### Create metadata table ####

sample_names <- colnames(raw_counts)
sampletype <- factor(c(rep("undiff_basal",4),rep("diff_basal",4),rep("diff_aicar",3), rep("diff_palm", 4)))
diff_state <- factor(c(rep("undiff", 4), rep("diff", 11)))
condition <- factor(c(rep("basal", 8), rep("aicar", 3), rep("palm", 4)))

metadata <- data.frame(sampleName = sample_names, group = sampletype, diff_state = diff_state, condition = condition)

# Ensure correct ordering
metadata$group <- relevel(metadata$group, ref = "undiff_basal")
metadata$diff_state <- relevel(metadata$diff_state, ref = "undiff")
metadata$condition <- relevel(metadata$condition, ref = "basal")

#### Create DESeq objects ####

# Validate alignment of counts and column data
if (!all(colnames(raw_counts) %in% metadata$sampleName)) stop("Sample names in countData and colData do not match.")

# Full DESeq dataset
dds <- DESeqDataSetFromMatrix(countData = raw_counts, colData = metadata, design = ~ group)
dds

# basal undiff vs. diff subset

metadata_basal <- metadata %>% filter(condition == "basal")
basal_counts <- raw_counts[, metadata_basal$sampleName]
dds_basal <- DESeqDataSetFromMatrix(countData = basal_counts, colData = metadata_basal, design = ~diff_state)

# differentiated conditions subset

metadata_diff <- metadata %>% filter(diff_state == "diff")
diff_counts <- raw_counts[, metadata_diff$sampleName]
dds_diff <- DESeqDataSetFromMatrix(countData = diff_counts, colData = metadata_diff, design = ~condition)

#### DESeq Analysis ####

# Pre-filtering
smallestGroupSize <- 4
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
