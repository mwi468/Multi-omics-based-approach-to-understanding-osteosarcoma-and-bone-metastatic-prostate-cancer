## plotting rfe model metrics
rm(list=ls())
library(caret)
library(dplyr)

#load("C:/Users/u/Documents/rfe_results.RData")

rf_stratified <- readRDS("rf_stratified_model.rds")

Combined_Annotations  <- read.delim("C:/Users/u/Dropbox/SharedDesktopFiles/Datasets/ProstateCancer_Human/CombinedAnnotations.tsv", header =TRUE, check.names = FALSE)
numeric_data_matrix <- read.delim("C:/Users/u/Dropbox/SharedDesktopFiles/Datasets/ProstateCancer_Human/CombinedNumerical_data.tsv", header =TRUE, check.names = FALSE)

unwanted_sites <- c("Other", "unknown", "Primary", "Prostate")

filtered_data <- Combined_Annotations %>%
  filter(!TISSUE_SITE %in% unwanted_sites)

## filter the numerical matrix based on this: 
sample_ids_to_keep <- filtered_data$SAMPLE_ID

# Filter columns in `numeric_data_matrix` to include only the sample IDs in `sample_ids_to_keep`

row_order <- match(sample_ids_to_keep, rownames(numeric_data_matrix))

# Filter the numeric_data_matrix to keep only the columns in sample_ids_to_keep and in the correct order
filtered_matrix <- numeric_data_matrix[row_order,]

filtered_data$LABEL <- ifelse(filtered_data$TISSUE_SITE == "Bone", 1, 0)
labels <- as.factor(filtered_data$LABEL)

# Perform RFE using your existing Random Forest model
rfe_control <- rfeControl(
  functions = rfFuncs,           # Use Random Forest functions
  method = "repeatedcv",         # Use repeated cross-validation
  number = 10,                   # Number of folds (matches stratified_control)
  repeats = 3,                   # Number of repeats (matches stratified_control)
  verbose = TRUE,                # Print progress
  saveDetails = TRUE,            # Save details for further analysis
  returnResamp = "final"         # Return final resampled results
)

# Perform RFE using the filtered_matrix and labels
rfe_results <- rfe(
  x = filtered_matrix,            # Feature matrix
  y = labels,                     # Class labels
  sizes = seq(10, 500, by = 10),  # Test subsets of features
  rfeControl = rfe_control,
  metric = "ROC",                 # Optimize for ROC (like stratified_control)
  method = "rf",                  # Random Forest
  ntree = rf_stratified$finalModel$ntree  # Use same number of trees as your existing RF model
)

# View RFE results
print(rfe_results)

save(rfe_results, file = "rfe_results.RData")
write.csv(rfe_data, "rfe_metrics.csv", row.names = FALSE)

rfe_data <- data.frame(
  num_features = rfe_results$results$Variables,
  accuracy = rfe_results$results$Accuracy,
  kappa = rfe_results$results$Kappa
)

par(mfrow = c(2, 1)) # 2 rows, 1 column layout


############# Plot rfe metrics 

# Plot full range
plot(rfe_data$num_features, rfe_data$accuracy, type = "b", pch = 19, col = "blue",
     xlab = "Number of Features", ylab = "Performance Metrics",
     main = "RFE Performance Metrics (Full Range)", ylim = c(0.85, 1))
lines(rfe_data$num_features, rfe_data$kappa, type = "b", pch = 17, col = "red")
legend("bottomright", legend = c("Accuracy", "Kappa"),
       col = c("blue", "red"), pch = c(19, 17), lty = 1)

rfe_data_zoom <- subset(rfe_data, num_features < 2000)

# Plot zoomed-in data (features below 1000)
plot(rfe_data_zoom$num_features, rfe_data_zoom$accuracy, type = "b", pch = 19, col = "blue",
     xlab = "Number of Features (Zoomed < 1000)", ylab = "Performance Metrics",
     main = "RFE Performance Metrics (Zoomed)", ylim = c(0.85, 1))
lines(rfe_data_zoom$num_features, rfe_data_zoom$kappa, type = "b", pch = 17, col = "red")
legend("bottomright", legend = c("Accuracy", "Kappa"),
       col = c("blue", "red"), pch = c(19, 17), lty = 1)


############# extract the top genes 

ranked_features <- rfe_results$variables

rfe_results$fit$importance

gene_importance <- rfe_results$fit$importance

# Sort the data based on the column for "1" (bone) in descending order
sorted_importance <- gene_importance[order(-gene_importance[, "1"]), ]

# View the top rows of the sorted importance scores
head(sorted_importance)
sorted_importance <- data.frame(sorted_importance)

# Get the unique feature names ranked by their importance
ranked_features <- ranked_features[order(ranked_features$Overall, decreasing = TRUE), ]

top_features <- (rownames(sorted_importance))

# Save the top 1000 features to a CSV file
write.csv(top_features, "topfeatures_rfe.csv", row.names = FALSE)

# Print the top 1000 features
cat("Top features selected by RFE:\n")
print(top_features)

genes_of_interest <- c("CX3CR1", "IL1B")

# Get the rank of these genes in the full ranked feature set
gene_ranks <- ranked_features[ranked_features$var %in% genes_of_interest, ]

# If ranked_features includes duplicate rankings across folds, aggregate to find average rank
gene_ranks <- gene_ranks %>%
  group_by(var) %>%
  summarize(Average_Rank = mean(Overall), .groups = "drop") %>%
  arrange(desc(Average_Rank))

# Display the ranks of the genes of interest
cat("Ranks of CX3CR1 and IL1B:\n")
print(gene_ranks)

# Check if they are in the top 1000
genes_in_top <- intersect(genes_of_interest, top_features)
if (length(genes_in_top) > 0) {
  cat("The following genes are in the top 1000 features:\n")
  print(genes_in_top)
} else {
  cat("None of the genes are in the top selected features.\n")
}


##################################### Plot heatmap in bone vs non bone using 200 top rfe genes

######## 1. Select  genes based on RFE
# 1. Extract top 200 RFE-selected genes
# 1. RFE-selected genes

top_200_genes <- top_features

# 2. Filter the gene expression matrix for the top 200 genes
top_genes_data <- filtered_matrix[, colnames(filtered_matrix) %in% top_200_genes]
top_genes_data <- top_genes_data[, top_200_genes]
rownames(top_genes_data) <- rownames(filtered_matrix)

# 3. Log-transform the gene expression data
#log_transformed_genes <- log2(top_genes_data + 1)

# 4. Standardize the data (Z-score normalization)
scaled_top_genes <- scale(top_genes_data, center = TRUE, scale = TRUE)

# 5. Reorder rows (samples) by Bone and Non-Bone labels
# Ensure Combined_Annotations matches the rownames of filtered_matrix
rownames(Combined_Annotations) <- Combined_Annotations$SAMPLE_ID  
annotation <- Combined_Annotations[rownames(scaled_top_genes), , drop = FALSE]  # Subset to match scaled_top_genes rownames

#Standardizing Tissue Column
annotation$Label <- ifelse(annotation$TISSUE_SITE == "Bone", "Bone", "Non-Bone")
annotation$TISSUE_SITE <- gsub("^(LN|Lymph node|Lymph Node)$", "Lymph Node", annotation$TISSUE_SITE, ignore.case = TRUE)  
annotation$TISSUE_SITE <- gsub("^(Soft tissue|Other Soft tissue)$", "Soft Tissue", annotation$TISSUE_SITE, ignore.case = TRUE)  

# Reorder samples to separate Bone and Non-Bone
bone_samples <- rownames(annotation[annotation$Label == "Bone", , drop = FALSE])
non_bone_samples <- rownames(annotation[annotation$Label == "Non-Bone", , drop = FALSE])
ordered_samples <- c(bone_samples, non_bone_samples)

# Reorder the matrix rows and annotation data
scaled_top_genes <- scaled_top_genes[ordered_samples, ]  # Reorder rows by sample order
annotation <- annotation[ordered_samples, , drop = FALSE]  # Reorder annotation data

# 6. Cap values in the range [-1, 1] for heatmap clarity
scaled_top_genes[scaled_top_genes > 1] <- 1
scaled_top_genes[scaled_top_genes < -1] <- -1

# 7. Define annotation colors for multiple annotations
annotation_colors <- list(
  Label = c("Bone" = "red", "Non-Bone" = "blue"),
  TISSUE_SITE = c("Bone" = "red","Lymph Node"= "green", "Liver" = "darkblue", "Soft Tissue" = "purple","Lung"="yellow","Adrenal"="brown","Skull base"="orange","Pelvic mass"="magenta","Brain"="cyan", "Primary Site"="lightyellow"),  # Replace with actual TISSUE_SITE values
  Subtype = c("NEPC" = "cyan", "Non_NEPC" = "brown", "unknown" = "gray"),  # Replace with actual Subtype values
  Study = c("SUC" = "blue", "UCSF" = "lightgreen", "WCM" = "orange")  # Replace with actual Study values
)

annotation_final <- annotation[, !(colnames(annotation) %in% c("SAMPLE_ID", "Tumor"))]

# 8. Generate the heatmap
library(pheatmap)
output_file <- "C:/Users/u/Documents/Capped_Z_Scored_Heatmap5.pdf"

pdf(output_file, width = 20, height = 20)  # Adjust width and height as needed

RdYlBu_colors <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(101)

pheatmap(
  mat = t(scaled_top_genes),  # Transpose for heatmap (genes as rows, samples as columns)
  annotation_col = annotation_final,  # Include additional annotations
  cluster_rows = FALSE,         # Do not cluster rows (keep RFE order)
  cluster_cols = FALSE,         # Do not cluster columns (manual order)
  show_rownames = TRUE,         # Show gene names
  show_colnames = FALSE,        # Hide sample names
  color = RdYlBu_colors,  # Gradient from blue to red
  breaks = seq(-1, 1, length.out = 101),  # Scale from -1 to 1
  annotation_colors = annotation_colors,  # Custom colors for annotations
  fontsize_row = 8,              # Smaller font size for row labels
  main = "Capped Z-Scored Heatmap of Top 200 RFE-Ranked Genes (Bone vs. Non-Bone)"
)

dev.off()

######### make tSNE using 750 top genes 

if (!require("Rtsne")) install.packages("Rtsne")
library(Rtsne)

## 750 with 14 works well

# Extract top 750 genes
top_200_genes <- unique(ranked_features$var)[1:140]

# Filter gene expression matrix
top_genes_data <- filtered_matrix[, colnames(filtered_matrix) %in% top_200_genes]
top_genes_data <- top_genes_data[, top_200_genes]
rownames(top_genes_data) <- rownames(filtered_matrix)

# Log-transform the gene expression data
#log_transformed_genes <- log2(top_genes_data + 1)

# Standardize the data (Z-score normalization)
scaled_top_genes <- scale(log_transformed_genes, center = TRUE, scale = TRUE)

# Remove WCM samples
filteredC <- Combined_Annotations[Combined_Annotations$Study != "WCM", ]

top_genes_data <- scaled_top_genes[rownames(scaled_top_genes) %in% filteredC$SAMPLE_ID, ]

filteredC <- filteredC[match(rownames(top_genes_data), filteredC$SAMPLE_ID), ]

# Check for NA values
if (anyNA(top_genes_data)) {
  stop("NA values detected in the dataset. Please handle missing values.")
}
filteredC$TISSUE_SITE <- gsub("^(LN|Lymph node|Lymph Node)$", "Lymph Node", filteredC$TISSUE_SITE, ignore.case = TRUE)  
filteredC$TISSUE_SITE <- gsub("^(Soft tissue|Other Soft tissue)$", "Soft Tissue",filteredC$TISSUE_SITE, ignore.case = TRUE)  

# t-SNE Analysis
tsne_result <- Rtsne(top_genes_data, dims = 2, perplexity = 12, max_iter = 1000)

# Plotting t-SNE results

output_file <- "C:/Users/u/Documents/tSNEScored_Heatmap3.pdf"
pdf(output_file, width = 10, height = 10)  # Adjust width and height as needed
par(mar = c(5, 4, 4, 10) + 0.1)  # Increase the right margin
# Plot the t-SNE result
plot(tsne_result$Y, col = as.factor(filteredC$TISSUE_SITE), pch = 19, 
     main = "t-SNE Visualization", xlab = "t-SNE 1", ylab = "t-SNE 2")

# Add the legend outside the plotting area on the right
legend("topright", inset = c(-.2, 0),  # Position legend outside the right margin
       legend = unique(filteredC$TISSUE_SITE), col = unique(as.factor(filteredC$TISSUE_SITE)), 
       pch = 19, xpd = TRUE)  # xpd 
dev.off()

### output summarry of top 140 genes
library(tible)

# Filter out WCM samples and update SAMPLE_ID
# Filter out WCM samples
filtered_data2 <- filtered_data %>%
  filter(Study != "WCM")

# Subset the matrix to match filtered_data2
filtered_matrix_selected <- filtered_matrix[rownames(filtered_matrix) %in% filtered_data2$SAMPLE_ID, ]

# Update SAMPLE_ID in filtered_data2 to match filtered_matrix_selected
filtered_data2 <- filtered_data2 %>%
  filter(SAMPLE_ID %in% rownames(filtered_matrix_selected))

# Proceed with gene summary creation
gene_summary <- filtered_matrix_selected %>%
  as.data.frame() %>%
  rownames_to_column(var = "SAMPLE_ID") %>%
  left_join(filtered_data2 %>% select(SAMPLE_ID, TISSUE_SITE), by = "SAMPLE_ID") %>%
  gather(key = "Gene", value = "Expression", -SAMPLE_ID, -TISSUE_SITE) %>%
  group_by(Gene) %>%
  summarize(
    Tissue_Site_Score = sum(TISSUE_SITE == "Bone"),
    Bone_Avg_Exp = mean(Expression[TISSUE_SITE == "Bone"], na.rm = TRUE),
    Non_Bone_Avg_Exp = mean(Expression[TISSUE_SITE != "Bone"], na.rm = TRUE),
    Up_or_Down_in_Bone = ifelse(
      mean(Expression[TISSUE_SITE == "Bone"], na.rm = TRUE) > mean(Expression[TISSUE_SITE != "Bone"], na.rm = TRUE),
      "Up",
      "Down"
    )
  ) %>%
  left_join(gene_scores, by = "Gene") %>%
  arrange(desc(Score))

gene_summary <- gene_summary %>%
  mutate(
    Bone_Avg_Exp = format(round(Bone_Avg_Exp, 2), nsmall = 2),
    Non_Bone_Avg_Exp = format(round(Non_Bone_Avg_Exp, 2), nsmall = 2),
    Score = format(round(Score, 4), nsmall = 4)
  )

# View the result
print(gene_summary)

write.csv(gene_summary, "gene_summaryUpdated.csv", row.names = FALSE)

