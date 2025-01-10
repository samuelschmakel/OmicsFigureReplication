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
lfcThreshold <- 0.6

res_diff_vs_undiff <- results(dds, contrast = c("condition", "diff", "undiff"), alpha = alpha, lfcThreshold = lfcThreshold)
res_AICAR_vs_undiff <- results(dds, contrast = c("condition", "AICAR", "undiff"), alpha = alpha, lfcThreshold = lfcThreshold)
res_palmitate_vs_undiff <- results(dds, contrast = c("condition", "palmitate", "undiff"), alpha = alpha, lfcThreshold = lfcThreshold)
res_diff_vs_AICAR <- results(dds, contrast = c("condition", "diff", "AICAR"), alpha = alpha, lfcThreshold = lfcThreshold)
res_diff_vs_palmitate <- results(dds, contrast = c("condition", "diff", "palmitate"), alpha = alpha, lfcThreshold = lfcThreshold)
res_AICAR_vs_palmitate <- results(dds, contrast = c("condition", "AICAR", "palmitate"), alpha = alpha, lfcThreshold = lfcThreshold)

#### PCA Plot ####
normalized_counts <- PCA_plotting(dds)

#### DEG list output ####
diff_vs_undiff_DEGs <- as.data.frame(res_diff_vs_undiff[res_diff_vs_undiff$padj < 0.05 & abs(res_diff_vs_undiff$log2FoldChange) > 0.6,])

diff_vs_undiff_DEG_counts <- as.data.frame(normalized_counts[rownames(normalized_counts) %in% rownames(diff_vs_undiff_DEGs),])

# Only keep samples that are involved in this comparison
diff_vs_undiff_DEG_counts <- diff_vs_undiff_DEG_counts[,1:8]

# Add direction of regulation
diff_vs_undiff_DEG_counts$regulation <- ifelse(diff_vs_undiff_DEGs$log2FoldChange > 0, "Upregulated", "Downregulated")

# Split matrix by regulation direction
upregulated <- rownames(diff_vs_undiff_DEG_counts[diff_vs_undiff_DEG_counts$regulation == "Upregulated", ])
downregulated <- rownames(diff_vs_undiff_DEG_counts[diff_vs_undiff_DEG_counts$regulation == "Downregulated", ])

# Reorder matrix
diff_vs_undiff_DEG_counts <- diff_vs_undiff_DEG_counts[c(upregulated, downregulated), ]

# Annotation for gene regulation
annotation_row <- data.frame(Regulation = diff_vs_undiff_DEG_counts$regulation)
rownames(annotation_row) <- rownames(diff_vs_undiff_DEG_counts)

# After ordering the data, remove the 'regulation' column before passing the matrix to pheatmap
diff_vs_undiff_DEG_counts_numeric <- diff_vs_undiff_DEG_counts[, -ncol(diff_vs_undiff_DEG_counts)]

# Create heatmap
pheatmap(
  diff_vs_undiff_DEG_counts_numeric, 
  cluster_rows = TRUE, 
  cluster_cols = TRUE, 
  annotation_row = annotation_row,
  scale = "row",  # Normalize rows to Z-scores
  color = colorRampPalette(c("blue", "white", "red"))(50),
  main = "Heatmap of DEGs"
)

#### GO Enrichment Analysis ####
upregulated <- rownames(diff_vs_undiff_DEG_counts[diff_vs_undiff_DEG_counts$regulation == "Upregulated",])
downregulated <- rownames(diff_vs_undiff_DEG_counts[diff_vs_undiff_DEG_counts$regulation == "Downregulated",])

# Create Entrez DEG List

# Remove version numbers from Ensembl Ids
up_genes <- gsub("\\..*$", "", upregulated)
down_genes <- gsub("\\..*$", "", downregulated)

# upregulated genes
entrez_ids <- bitr(up_genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

entrez_gene_list <- entrez_ids$ENTREZID

go_results <- enrichGO(
  gene          = entrez_gene_list, 
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",             # Ontology: "BP" (Biological Process), "MF" (Molecular Function), or "CC" (Cellular Component)
  pAdjustMethod = "BH",             # Adjust p-values for multiple testing
  pvalueCutoff  = 0.05, 
  qvalueCutoff  = 0.05,
  readable      = TRUE              # Convert Entrez IDs to gene symbols
)

# Visualize results
barplot(go_results, showCategory = 10)  # Top 10 enriched categories

dotplot(go_results, showCategory = 10)
cnetplot(go_results, showCategory = 5)

# downregulated genes
entrez_ids <- bitr(down_genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

entrez_gene_list <- entrez_ids$ENTREZID

go_results <- enrichGO(
  gene          = entrez_gene_list, 
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",             # Ontology: "BP" (Biological Process), "MF" (Molecular Function), or "CC" (Cellular Component)
  pAdjustMethod = "BH",             # Adjust p-values for multiple testing
  pvalueCutoff  = 0.05, 
  qvalueCutoff  = 0.05,
  readable      = TRUE              # Convert Entrez IDs to gene symbols
)

# Visualize results
barplot(go_results, showCategory = 10)  # Top 10 enriched categories

dotplot(go_results, showCategory = 10)
cnetplot(go_results, showCategory = 5)