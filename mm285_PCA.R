PCA analysis of mm285 syntenic probes 
library(tidyverse)
library(readxl)
library(Rtsne)
library(factoextra)
library(pheatmap)
library(gridExtra)

# Create a data frame for plotting
plot_data <- data.frame(
  PC1 = pca_res$x[, p1],
  PC_Avg_2_10 = rowMeans(pca_res$x[, 2:10], na.rm = TRUE),  # Average of PCs 2-10
  TISSUE = Exp$TISSUE[match(rownames(betas), Exp$SampleID)],
  SPECIES = Exp$SPECIES[match(rownames(betas), Exp$SampleID)]
)

# Basic PCA plot
ggplot(plot_data, aes(x = PC1, y = PC_Avg_2_10, color = TISSUE, shape = SPECIES)) +
  geom_point(size = 2) +
  labs(title = "PCA of Methylation Data", x = paste("PC", p1), y = "Average of PCs 2-10") +
  theme_minimal() +
  theme(legend.position = "right")


# Getting variance eigenvalues for PCs used in plot
x <- get_eig(pca_res)

# Principal Component Correlation with Sample Types
df <- pca_res[["x"]]

# Species difference analysis
test <- subset(df, rownames(df) %in% Exp[Exp$SPECIES == 'Human',]$SampleID)
test2 <- subset(df, rownames(df) %in% Exp[Exp$SPECIES != 'Human',]$SampleID)
edf <- data.frame(matrix(NA, nrow = 9, ncol = ncol(test)))

for (i in 1:ncol(test)) {
  rownames(edf)[1] <- 'Species'
  colnames(edf)[i] <- colnames(test)[i]
  p <- wilcox.test(test[, i], test2[, i])
  edf[1, i] <- log(p[["p.value"]], 10)
}

# Tissue difference analysis
j <- unique(Exp$TISSUE)

for (i in 1:length(j)) {
  test <- subset(df, rownames(df) %in% Exp[Exp$TISSUE == j[i],]$SampleID)
  test2 <- subset(df, rownames(df) %in% Exp[Exp$SPECIES != j[i],]$SampleID)
  rownames(edf)[1 + i] <- j[i]

  for (z in 1:ncol(test)) {
    p <- wilcox.test(test[, z], test2[, z])
    edf[i + 1, z] <- log(p[["p.value"]], 10)
  }
}

# Only keep the first 10 Principal components
final_edf <- edf[, 1:10]

# Scree plot
scree_plot <- fviz_eig(pca_res)

# Heatmap
heatmap_plot <- pheatmap(final_edf, cluster_cols = FALSE, cluster_rows = FALSE,
                          cellheight = 15, cellwidth = 30, main = "PCA Heatmap",
                          show_colnames = TRUE, show_rownames = TRUE, fontsize = 10,
                          legend = TRUE, legend_labels = c('HighPv', 'UpRegulated'),
                          height = 5, width = 10)

# Combine and save plots
combined_plot <- grid.arrange(arrangeGrob(scree_plot, heatmap_plot[[4]], nrow = 2, ncol = 1))
ggsave("~/gallery/20210811_SyntenicProbes_PCAheatmap_LOG10.pdf", combined_plot, width = 5, height = 5, dpi = 300, scale = 1.5)


