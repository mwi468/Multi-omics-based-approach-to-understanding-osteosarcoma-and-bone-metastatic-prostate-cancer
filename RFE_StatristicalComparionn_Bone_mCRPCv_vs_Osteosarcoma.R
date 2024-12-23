rm(list=ls())

library(dplyr)
library(reshape2)
library(hgnc)
library(rtracklayer)

###########################################3
####### Functions for proper formatting
aggregate_transcript_to_gene_level <- function(data) {
  # Ensure there is at least one gene column and one sample column
  if (ncol(data) < 2) {
    stop("Data must have at least one gene column and one sample column with numeric values.")
  }
  
  # Rename the first column to "gene_id" for consistency
  colnames(data)[1] <- "gene_id"
  
  # Aggregate values by summing across numeric entries for each gene
  gene_level_data <- data %>%
    group_by(gene_id) %>%
    summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE))) %>%
    ungroup()
  
  return(gene_level_data)
}

##################################
####### Get RFE genes

genes <- read.csv("C:/Users/u/Documents/gene_summaryUpdated.csv")
#genes <- read.csv("C:/Users/Waleed/Dropbox/SharedDesktopFiles/gene_summaryUpdated.csv")

genes <- genes[trimws(genes$Score) != "NA", ]
genes <- genes[genes$Score !="NA" ,]

genes$rank <- rownames(genes)
  
## fix the one gene 
# Import the latest HGNC dataset
#hgnc_data <- import_hgnc_dataset()
#hgnc_data <- import_hgnc_dataset(file = "https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt")

# Create a named vector for mapping previous symbols to current symbols
#symbol_map <- setNames(hgnc_data$symbol, hgnc_data$prev_symbol)

# Update your gene symbols
#genes$genes <- ifelse(genes$genes %in% names(symbol_map),
#                      symbol_map[genes$genes],
#                      genes$genes)

##### Combining Primary Bone datasets
### Combine Primary Bone Datasets from TARGET and Shanghai Hospital

Bone1 <- read.table("C:/Users/u/Dropbox/SharedDesktopFiles/TARGET_geneExpressionMerged.tsv", header =TRUE, check.names = FALSE)
#Bone <- read.table("C:/Users/Waleed/Dropbox/SharedDesktopFiles/TEST.tsv", header =TRUE)

samples <- read.delim("C:/Users/u/Dropbox/TARGET/TARGET_GeneExpressionSample sheet.tsv")

file_names <- sub("\\.rna_seq\\.augmented_star_gene_counts\\.tsv$", "", samples$File.Name)  # Remove extensions
sample_ids <- samples$Case.ID  # Corresponding sample IDs

# Create a named vector for mapping: file name -> sample ID
name_map <- setNames(sample_ids, file_names)

# Get the column names of Bone
bone1_colnames <- colnames(Bone1)

# Replace column names in Bone using the mapping
updated_colnames <- sapply(bone1_colnames, function(col) {
  if (col %in% names(name_map)) {
    name_map[col]  # Replace with corresponding sample ID
  } else {
    col  # Keep original name if no match
  }
})

# Update the column names in Bone
colnames(Bone1) <- updated_colnames

Bone1 <- Bone1[,2:ncol(Bone1)]
names(Bone1)[1] <- 'Geneid'

# Verify the updated column names
head(colnames(Bone1))

Bone2 <- read.delim("e:/ShanghaiOsteosarcoma_WES_RNAseq/TPM/TPMsorted/final_combined.txt")

## convert ensmbl ids to gene names
gtf_file <- "C:/Users/u/Dropbox/SharedDesktopFiles/gencode.v47.annotation.gtf"
# Import the GTF file
gtf_data <- rtracklayer::import(gtf_file)
# Convert to a data frame for easy handling
gtf_df <- as.data.frame(gtf_data)

# Preview the data
head(gtf_df)

mapping <- gtf_df[, c("gene_id", "gene_name")]

# Remove duplicate mappings
mapping <- unique(mapping)

# Preview the mapping
head(mapping)

colnames(mapping) <- c("Geneid", "Gene_Name")
merged_data <- merge(Bone2, mapping, by = "Geneid", all.x = TRUE)

# Replace Ensembl IDs with gene names if available
merged_data$Geneid <- ifelse(!is.na(merged_data$Gene_Name), merged_data$Gene_Name, merged_data$Geneid)
aggregated_Bone1_data <- aggregate_transcript_to_gene_level(Bone1)
aggregated_Bone2_data <- aggregate_transcript_to_gene_level(merged_data)
common_genes <- intersect(aggregated_Bone1_data$gene_id, aggregated_Bone2_data$gene_id)
length(common_genes)
subset1 <-aggregated_Bone1_data[aggregated_Bone1_data$gene_id %in% common_genes, ]
nrow(subset1)
subset2 <- aggregated_Bone2_data[aggregated_Bone2_data$gene_id %in% common_genes, ]
nrow(subset2)

combined_data <- merge(subset1, subset2, by = "gene_id", all.x = TRUE)
nrow(combined_data)
length(unique(combined_data$Geneid))

saveRDS(combined_data, "C:/Users/u/Dropbox/SharedDesktopFiles/Bone_CombinedGeneExpression.rds")

## subset for only RFE genes 

rownames(combined_data) <- combined_data$gene_id
combined_data <- t(combined_data[,2:ncol(combined_data)])
combined_data <- data.frame(combined_data, check.names = FALSE)
combined_data <- combined_data[, colnames(combined_data) %in% genes$Gene]
saveRDS(combined_data, "C:/Users/u/Dropbox/SharedDesktopFiles/Bone_RFEsubet.rds")

######################
# Format the Combined Prostate Cancer and Combined Primary Osteosarcoma 

## load prostate cancer dataset
numeric_data_matrix <- read.delim("C:/Users/u/Dropbox/SharedDesktopFiles/Datasets/ProstateCancer_Human/CombinedNumerical_data.tsv", header =TRUE, check.names = FALSE)
#numeric_data_matrix <- read.delim("C:/Users/Waleed/Dropbox/SharedDesktopFiles/Datasets/ProstateCancer_Human/CombinedNumerical_data.tsv", header =TRUE, check.names = FALSE)

Combined_Annotations  <- read.delim("C:/Users/u/Dropbox/SharedDesktopFiles/Datasets/ProstateCancer_Human/CombinedAnnotations.tsv", header =TRUE, check.names = FALSE)
#Combined_Annotations  <- read.delim("C:/Users/Waleed/Dropbox/SharedDesktopFiles/Datasets/ProstateCancer_Human/CombinedAnnotations.tsv", header =TRUE, check.names = FALSE)

filtered <- Combined_Annotations[Combined_Annotations$Study != "WCM",]
filtered <- filtered[filtered$TISSUE_SITE == 'Bone',]

ProstateNumericMatrix <- numeric_data_matrix[rownames(numeric_data_matrix) %in% filtered$SAMPLE_ID, ]
ProstateNumericMatrix <- ProstateNumericMatrix[, colnames(ProstateNumericMatrix)  %in% colnames(combined_data)]

common_genes <- intersect(colnames(ProstateNumericMatrix), colnames(combined_data))

# Subset and reorder columns in both data frames based on common genes
ProstateAligned <- ProstateNumericMatrix[, common_genes]
BoneAligned <- combined_data[, common_genes]

#### Identifying genes with similar expression
combined_data_all <- rbind(ProstateAligned,BoneAligned)

# Ensure group labels match rows in combined_data
group_labels <- c(rep("Prostate", nrow(ProstateAligned)), rep("Bone", nrow(BoneAligned)))
group_labels <- factor(group_labels, levels = c("Prostate", "Bone"))

# Confirm the levels
levels(group_labels)  # Should return "Prostate" "Bone"

# Create a design matrix for the analysis
design <- model.matrix(~ group_labels)

library(limma)

# Create a design matrix for the analysis
design <- model.matrix(~ group_labels)

# Fit a linear model to the data
fit <- lmFit(t(combined_data_all), design)

# Apply empirical Bayes moderation
fit <- eBayes(fit)

# Extract the results for the comparison between groups
results <- topTable(fit, coef = 2, number = Inf, adjust.method = "BH")

# Bone is reference group
# Pos. FC: Indicates upregulation in the "Bone" group compared to the "Prostate" group.
# Neg. FC downregulation in Bone group 

results$PrimaryBone <- ifelse(results$logFC > 0, "Up", "Down")
########## Similar expression: Limma did not find sig diff expression
p_value_threshold <- 0.05
# Identify genes with adjusted p-values greater than the threshold
similar_genes <- results[results$adj.P.Val > p_value_threshold, ]
similar_genes$PrimaryBone <- 'Similar'

## Down regulated genes 
## Genes found to be downregulated in prostate bone metastases is even more downregulated

## down genes 
genes_down <- genes[genes$Up_or_Down_in_Bone == "Down",]
genes_down <- genes_down$Gene
Down_genes <- results[rownames(results) %in% genes_down,]
Down_genes <- Down_genes[Down_genes$PrimaryBone == "Down",]

## Up genes 
genes_up <- genes[genes$Up_or_Down_in_Bone == "Up",]
genes_up <- genes_up$Gene
Up_genes <- results[rownames(results) %in% genes_up,]
Up_genes <- Up_genes[Up_genes$PrimaryBone == "Up",]

## add RFE ranking
CombinedGenes <- rbind(Up_genes,Down_genes)
## there might be up or down according to FC but not sig accoring to limma remove these before adding
## otherwise multiple rows [bc we marked them similar and the y will be marked up or down]
CombinedGenes <- CombinedGenes[!rownames(CombinedGenes) %in% rownames(similar_genes),]
CombinedGenes <- rbind(CombinedGenes,similar_genes)
CombinedGenes$Gene <- rownames(CombinedGenes)

genes <- genes %>%
  rename(MestastaticBone = Up_or_Down_in_Bone)

merged_df <- merge(CombinedGenes, genes[, c("Gene","MestastaticBone","rank")], by = "Gene", all.x = TRUE)

## visualizing table subsets
## upregulated genes 
sub <- merged_df[merged_df$MestastaticBone == "Up",]
sub <- sub[order(sub$PrimaryBone, decreasing=FALSE), ]

#cat(paste(round(sub$logFC, digits =2), collapse = "\n"))
cat(paste(sub$Gene, collapse = "\n"))
#cat(paste(round(sub$adj.P, digits =15), collapse = "\n"))
cat(paste(format(sub$adj.P.Val, scientific = TRUE, digits = 3), collapse = "\n"))

saveRDS(merged_df , "C:/Users/u/Dropbox/SharedDesktopFiles/RFE_geneComparisonTable.rds")
