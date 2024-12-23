rm(list=ls())

# Load necessary libraries
library(ggplot2)
library(matrixStats)
library(data.table)
library(dplyr)

# Load the data
x <- read.table("C:/Users/Waleed/Downloads/site_wcdt.csv", sep =",", header = TRUE)
xa <- read.table("C:/Users/Waleed/Dropbox/SharedDesktopFiles/Lundberg_at_al/counts_and_TPMs/Counts_TPM_Genes.txt", header = TRUE, check.names = FALSE)
xb <- read.table("C:/Users/Waleed/Dropbox/SharedDesktopFiles/Lundberg_at_al/counts_and_TPMs/Counts_TPM_Genes_PROMOTE.txt", header = TRUE, check.names = FALSE)

# Merge the datasets by 'FEATURE_ID'
x1 <- merge(xa[1:(ncol(xa)-1)], xb[1:(ncol(xb))], by="FEATURE_ID")
rm(xa, xb)

# Define the genes of interest
Stemnessgenes <- c('SALL4', 'POU5F1', 'ALDH1A1', 'NR5A2', 'SOX2', 'NANOG', 'CX3CR1')
EMT_genes <- c('ZEB1', 'SNAI1', 'SNAI2', 'TWIST1', 'CDH2', 'IL1B', 'CX3CR1')

# Subset and transform the data for stemness genes
x2_stemness <- x1[x1$gene_name %in% Stemnessgenes, ]
t1_stemness <- data.frame(t(x2_stemness[, 2:(ncol(x2_stemness)-1)]))
colnames(t1_stemness) <- x2_stemness$gene_name

# Log2 transform the data with a pseudocount to handle zero values
n.exp_stemness <- log2(t1_stemness + 0.001)

# Standardize the data (Z-score normalization)
zexp_stemness <- data.frame(scale(n.exp_stemness, center = TRUE, scale = TRUE))

# Calculate the mean of the stemness genes
zexp_stemness$mean_stemness <- rowMeans(zexp_stemness[, Stemnessgenes], na.rm = TRUE)

# Subset and transform the data for EMT genes
x2_EMT <- x1[x1$gene_name %in% EMT_genes, ]
t1_EMT <- data.frame(t(x2_EMT[, 2:(ncol(x2_EMT)-1)]))
colnames(t1_EMT) <- x2_EMT$gene_name

# Log2 transform the data with a pseudocount to handle zero values
n.exp_EMT <- log2(t1_EMT + 0.001)

# Standardize the data (Z-score normalization)
zexp_EMT <- data.frame(scale(n.exp_EMT, center = TRUE, scale = TRUE))

# Calculate the average of EMT genes to represent EMT activity
zexp_EMT$EMT_Activity <- rowMeans(zexp_EMT[, EMT_genes], na.rm = TRUE)

# Load the immune-filtered dataset and subset rows
df <- readRDS("C:/Users/Waleed/Dropbox/SharedDesktopFiles/CX3CR1/ImmuneFiltered3.rds")
