# helper functions

get_results <- function(name, normalized_counts, alpha = 0.05, lfcThreshold = 1) {
  res <- results(dds, name=name, alpha=alpha, lfcThreshold=lfcThreshold)
  return (res)
}

PCA_plotting <- function(dds) {
  rld <- rlog(dds, blind = TRUE)
  
  normalized_counts <- assay(rld)
  pca <- prcomp(t(normalized_counts))
  percentVar <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 1)
  
  # Create a data frame for ggplot
  pca_df <- data.frame(
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    condition = colData(dds)$condition
  )
  
  # Figure 2a
  pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = condition, shape = condition)) +
    geom_point(size = 3, color = "black", aes(fill = condition)) +
    scale_shape_manual(values = c(24, 21, 21, 21)) +
    scale_fill_manual(values = c("red", "yellow", "green", "blue"))  +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    ggtitle("PCA of Gene Expression by Condition") +
    theme_bw()
  
  # Save the results into the results folder
  if (!dir.exists("results")) {
    dir.create("results")
  }
  
  ggsave(filename = "results/figure2a.png", plot = pca_plot, width = 8, height = 6, dpi = 300)
  
  return (normalized_counts)
}

get_counts <- function(res, normalized_counts, alpha, lfcThreshold) {
  DEGs <- as.data.frame(res[res$padj < alpha & abs(res$log2FoldChange) > lfcThreshold,])
  
  return (as.data.frame(normalized_counts[rownames(normalized_counts) %in% rownames(DEGs),]))
}

reorder_counts <- function(cts, DEG_info) {
  # Add direction of regulation
  DEG_info <- DEG_info[rownames(cts),]
  
  cts$regulation <- ifelse(DEG_info$log2FoldChange > 0, "Upregulated", "Downregulated")
  
  # Split matrix by regulation direction
  upregulated <- rownames(cts[cts$regulation == "Upregulated", ])
  downregulated <- rownames(cts[cts$regulation == "Downregulated", ])
  
  # Reorder matrix
  cts <- cts[c(upregulated, downregulated), ]
  
  # Annotation for gene regulation
  annotation_row <- data.frame(Regulation = cts$regulation)
  rownames(annotation_row) <- rownames(cts)
  
  # After ordering the data, remove the 'regulation' column before passing the matrix to pheatmap
  cts_numeric <- cts[, -ncol(cts)]
  
  return (list(cts = cts_numeric, annotation = annotation_row))
}

GO_enrichment <- function(counts, name) {
  upregulated <- rownames(counts[cts$regulation == "Upregulated",])
  downregulated <- rownames(counts[cts$regulation == "Downregulated",])
  
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
  
  bplot_up <- barplot(go_results, showCategory = 10)
  
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
  
  bplot_down <- barplot(go_results, showCategory = 10)
  
  #dotplot(go_results, showCategory = 10)
  #cnetplot(go_results, showCategory = 5)
  
  # ensure the directory exists, if it doesn't, create it
  if (!dir.exists("results")) {
    dir.create("results")
  }
  
  if (name == "diff_vs_undiff") {
    ggsave(filename = "results/figure2c_diff_GOterms.png", plot = bplot_up, width = 8, height = 6, dpi = 300)
    ggsave(filename = "results/figure2c_undiff_GOterms.png", plot = bplot_down, width = 8, height = 6, dpi = 300)
  } else if (name == "diff_vs_AICAR") {
    ggsave(filename = "results/figure2e_diff_GOterms.png", plot = bplot_up, width = 8, height = 6, dpi = 300)
    ggsave(filename = "results/figure2e_AICAR_GOterms.png", plot = bplot_down, width = 8, height = 6, dpi = 300)
  } else if (name == "diff_vs_palmitate") {
    ggsave(filename = "results/figure2f_diff_GOterms.png", plot = bplot_up, width = 8, height = 6, dpi = 300)
    ggsave(filename = "results/figure2f_palmitate_GOterms.png", plot = bplot_down, width = 8, height = 6, dpi = 300)
  } else {
    print("Invalid comparison name")
  }
}