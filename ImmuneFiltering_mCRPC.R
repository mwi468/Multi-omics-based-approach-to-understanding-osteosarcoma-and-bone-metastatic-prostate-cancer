rm(list=ls())

library(tidyr)
library(readxl)
library(pheatmap)
library(RColorBrewer)
library(matrixStats)
library(ggplot2)
library(dplyr)

# Read data
xa <- read.table("C:/Users/Waleed/Dropbox/SharedDesktopFiles/Lundberg_at_al/counts_and_TPMs/Counts_TPM_Genes.txt", header =TRUE, check.names = FALSE)
xb <- read.table("C:/Users/Waleed/Dropbox/SharedDesktopFiles/Lundberg_at_al/counts_and_TPMs/Counts_TPM_Genes_PROMOTE.txt", header =TRUE, check.names = FALSE)
exp <- merge(xa[1:(ncol(xa)-1)], xb[1:(ncol(xb))], by="FEATURE_ID")
rm(xa, xb)

# Read immune gene list
g <- read_excel("C:/Users/Waleed/Dropbox/SharedDesktopFiles/ImmuneGenes_PMID_29230012.xlsx")
g <- tidyr::separate_rows(g, Required_Genes, sep=",")
g <- g[order(g$Cell_type),]

# Filter expression data for required genes
exp2 <- exp[exp$gene_name %in% g$Required_Genes,]
reorder_idx <- match(g$Required_Genes, exp2$gene_name)
exp2 <- exp2[reorder_idx,]
rownames(exp2) <- exp2$gene_name
expF <- data.frame(t(exp2[, 2:(ncol(exp2)-1)]))

# Replace zeros and filter genes
expF_filtered <- expF[, colSums(expF == 0) < (0.5 * nrow(expF))]

# Log transformation
expF_filtered <- log2(expF_filtered + .001)

# Z-score calculation
standard_zscore <- function(x) {
  m <- mean(as.numeric(x))
  sd_value <- sd(as.numeric(x))
  if (sd_value == 0) sd_value <- 0.001
  n <- (x - m) / sd_value
  return(n)
}
zexp <- data.frame(apply(expF_filtered, 2, standard_zscore))
rownames(zexp) <- rownames(expF_filtered)

# Create annotation
ann <- data.frame(unique(names(zexp)))
names(ann)[1] <- names(g)[2]
ann <- merge(ann, g, by = names(g)[2], all.x = FALSE, no.dups = TRUE)
rownames(ann) <- ann$Required_Genes
ann <- subset(ann, select = "Cell_type")

# Plot heatmap
p <- pheatmap(zexp, colorRampPalette(rev(brewer.pal(n = 6, name = "RdBu")))(6), breaks = c(-3, -2, -1, 0, 1, 2, 3), cluster_rows = FALSE, cluster_cols = FALSE, cellheight = 2, cellwidth = 8, cex = .75, show_colnames = TRUE, show_rownames = FALSE, annotation_col = ann, treeheight_row = 0, treeheight_col = 0, fontsize = 7, legend = TRUE, height = 10, width = 5)
ggsave("C:/Users/Waleed/Dropbox/SharedDesktopFiles/ImmuneContamination_Martin4.pdf", plot = p, width = 10, height = 15, limitsize = FALSE)

# Filter by immune cell types and calculate averages
calculate_cell_type_average <- function(cell_type, zexp, g) {
  genes <- g[g$Cell_type == cell_type,]$Required_Genes
  if (length(genes) > 0) {
    cell_df <- data.frame(zexp[, colnames(zexp) %in% genes])
    cell_df$Average <- rowMedians(as.matrix(cell_df))
    return(cell_df$Average)
  }
  return(NULL)
}

# Immune cell types to include (excluding Endothelial Cell and Cancer Associated Fibroblast)
immune_cell_types <- c("Tcells", "CD4+ Tcells", "CD8a+ Tcells", "CD8b+ Tcells", "Tregs", "B cell", "Macrophage/Monocyte", "Dendritic Cell", "Natural Killer Cell")
NewDF <- sapply(immune_cell_types, function(cell_type) calculate_cell_type_average(cell_type, zexp, g))
NewDF <- data.frame(NewDF)
colnames(NewDF) <- immune_cell_types
rownames(NewDF) <- rownames(zexp)

# Plot heatmap for averages
p <- pheatmap(NewDF, colorRampPalette(rev(brewer.pal(n = 6, name = "RdBu")))(6), breaks = c(-3, -2, -1, 0, 1, 2, 3), cluster_rows = FALSE, cluster_cols = FALSE, cellheight = 2, cellwidth = 8, cex = 1, show_colnames = TRUE, show_rownames = FALSE, treeheight_row = 0, treeheight_col = 0, fontsize = 7, legend = TRUE, height = 10, width = 5)
#ggsave("C:/Users/Waleed/Dropbox/SharedDesktopFiles/ImmuneContaminationSimplified2.pdf", plot = p ,width = 10, height = 15)

# Filter by immune cell contamination levels
test <- rowSums(NewDF[,1:9] >= 1.25)
filtered_rows <- NewDF[test == 0, ]

# Remove samples with lesser contribution from multiple immune cells
num_cols <- ncol(filtered_rows) - 1  # Excluding the last column
threshold <- 0.75 * num_cols
rows_to_keep <- apply(filtered_rows[, 1:num_cols] > 0.5, 1, function(row) sum(row) <= threshold)
final_filtered_rows <- filtered_rows[rows_to_keep, ]

# Save final filtered data
saveRDS(final_filtered_rows, "C:/Users/Waleed/Dropbox/SharedDesktopFiles/CX3CR1/ImmuneFiltered3.rds")
