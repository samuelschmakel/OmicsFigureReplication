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

GO_enrichment <- function(counts, name) {
  upregulated <- rownames(counts[diff_vs_undiff_DEG_counts$regulation == "Upregulated",])
  downregulated <- rownames(counts[diff_vs_undiff_DEG_counts$regulation == "Downregulated",])
  
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