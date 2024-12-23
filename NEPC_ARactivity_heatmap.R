rm(list=ls())

# Load necessary libraries
library(data.table)
library(stringr)
library(RColorBrewer)
library(pheatmap)

# Load data
exp <- read.table("C:/Users/Waleed/Dropbox/SharedDesktopFiles/IL1B_collaboration_data/featurecounts_tpm_gene_gencodev28_n100.tsv", header = TRUE)
list <- read.table("C:/Users/Waleed/Dropbox/SharedDesktopFiles/NEPC genes.txt", sep = " ")
names <- fread("C:/Users/Waleed/Dropbox/SharedDesktopFiles/gencode.v28.simplified.txt", header = FALSE)

# AR activity genes 
argenes <- c('KLK3','KLK2', 'FKBP5', 'STEAP1', 'STEAP2', 'PLPP1', 'RAB3B', 'NKX3-1', 'ACSL3')

# Subset expression for relevant genes 
names2keep1 <- names[names$V3 %in% list$V1,]
names(names2keep1)[2] <- 'feature_id'
names2keep2 <- names[names$V3 %in% argenes,]
names(names2keep2)[2] <- 'feature_id'

# Subset exp based on ensembl id in names2keep
exp1 <- exp[exp$feature_id %in% names2keep1$feature_id,]
exp1 <- merge(names2keep1, exp, by = 'feature_id' )
exp2 <- exp[exp$feature_id %in% names2keep2$feature_id,]
exp2 <-  merge(names2keep2, exp, by = 'feature_id' )
exp <- rbind(exp1,exp2)

# Split column separating gene id and names
n.exp <- exp[,4:ncol(exp)]
rownames(n.exp) <- exp$V3

# Convert 0 to a small number to avoid issues in log transformations
n.exp[n.exp == 0] = .01

# Z-score normalization across rows (genes)
ZscoreRow <- function(row) {
  (row - mean(row, na.rm = TRUE)) / sd(row, na.rm = TRUE)
}
zexp <- t(apply(n.exp, 1, ZscoreRow))
zexp <- data.frame(zexp)
rownames(zexp) <- rownames(n.exp)
colnames(zexp) <- colnames(n.exp)

zexp <- data.frame(t(zexp))

# Calculate AR and NEPC activity as the average Z-score of respective gene sets

zexp$ARActivity <- rowMeans(zexp[,(which(colnames(zexp)=="KLK3"):which(colnames(zexp)=="ACSL3"))], na.rm = TRUE)
zexp$NEPCActivity <- rowMeans(zexp[,(which(colnames(zexp) %in% list$V1))], na.rm = TRUE)

# Reorder based on AR activity
zexp <- zexp[order(zexp$ARActivity),]

# Define color palette with custom thresholds
my_breaks <- c(-1.5, -1, 0, 1, 1.5)
my_colors <- colorRampPalette(c("blue", "white", "red"))(length(my_breaks) - 1)

# Plot heatmap with custom color breaks
pheatmap(zexp[,1:ncol(zexp)],
         scale = "none",
         color = my_colors,
         breaks = my_breaks,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         cellheight = 4,
         cellwidth = 8,
         cex = 1,
         show_colnames= TRUE,
         show_rownames= FALSE,
         treeheight_row = 0,
         treeheight_col = 0,
         fontsize = 7,
         legend = TRUE,
         height = 10,
         width = 15)

correlation_pearson <- cor.test(zexp$NEPCActivity, zexp$ARActivity, method = 'pearson')
correlation_spearman <- cor.test(zexp$NEPCActivity, zexp$ARActivity, method = 'spearman')

# Print correlation results
print(correlation_pearson)
print(correlation_spearman)
