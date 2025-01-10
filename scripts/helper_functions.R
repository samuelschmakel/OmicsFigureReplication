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