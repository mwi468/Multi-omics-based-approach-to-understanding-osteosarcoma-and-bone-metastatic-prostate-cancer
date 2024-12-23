rm(list=ls())

library(dplyr)
library(reshape2)
library(survival)
library(ggplot2)

## load the formatted bone gene expression data set 

bones <- readRDS("C:/Users/u/Dropbox/SharedDesktopFiles/Bone_CombinedGeneExpression.rds")
RFEgenes <- readRDS("C:/Users/u/Dropbox/SharedDesktopFiles/RFE_geneComparisonTable.rds")

#subset the bone gene expression for only the RFEmatching genes (from comparing primary bone and metastatic bone)

bones <- bones[bones$gene_id %in% RFEgenes$Gene,]

rownames(bones) <- bones$gene_id
bones <- bones[,2:ncol(bones)]
bones <- data.frame(t(bones))
############### Clinical

## load clinical data for OS

# Step 1: Load clinical data
clinical_file1 <- "C:/Users/u/Dropbox/TARGET/clinical.project-target-os.2024-11-25/clinical.tsv"
clinical_file1 <- read.delim(clinical_file1, header = TRUE, stringsAsFactors = FALSE)

# Subset relevant columns
clinical_data_main <- clinical_file1[, c("case_submitter_id", "days_to_death", "days_to_last_follow_up", "vital_status", "age_at_index", "days_to_birth", "age_at_diagnosis", "days_to_diagnosis","year_of_diagnosis")]
clinical_data_main <- clinical_data_main %>% distinct()

followUP <- "C:/Users/u/Dropbox/TARGET/clinical.project-target-os.2024-11-25/follow_up.tsv"
follow_up_data <- read.delim(followUP , header = TRUE, stringsAsFactors = FALSE)
follow_up_data_main <- follow_up_data[,  c("case_submitter_id", "days_to_first_event", "days_to_follow_up", "days_to_progression", "first_event","year_of_follow_up")]

## merge information for duplicates

follow_up_data_main <- follow_up_data_main %>%
  mutate(across(where(is.character), ~na_if(.x, "'--"))) %>%  # Replace "'--" in character columns
  mutate(across(where(is.numeric), ~replace(.x, .x == as.numeric("'--"), NA))) # Replace in numeric columns, if necessary

clinical_data_main <- clinical_data_main %>%
  mutate(across(where(is.character), ~na_if(.x, "'--"))) %>%  # Replace "'--" in character columns
  mutate(across(where(is.numeric), ~replace(.x, .x == as.numeric("'--"), NA))) # Replace in numeric columns, if necessary


follow_up_data_main <- follow_up_data_main %>%
  group_by(case_submitter_id) %>%
  summarise(
    days_to_first_event = first(na.omit(days_to_first_event)),
    days_to_follow_up = first(na.omit(days_to_follow_up)),
    days_to_progression = first(na.omit(days_to_progression)),
    first_event = first(na.omit(first_event)),
    year_of_follow_up = first(na.omit(year_of_follow_up))
  ) %>%
  ungroup()

merged_data <- clinical_data_main %>%
  left_join(follow_up_data_main, by = "case_submitter_id")
merged_data <- merged_data[merged_data$case_submitter_id %in% rownames(bones),]

merged_data_data_main <- merged_data %>%
  mutate(
    # Derive survival_time based on vital_status
    survival_time = case_when(
      vital_status == "Dead" ~ as.numeric(days_to_death),
      vital_status == "Alive" ~ as.numeric(days_to_follow_up),
      TRUE ~ NA_real_  # For "Not Reported" or other cases
    ),
    # Derive event status
    event = case_when(
      vital_status == "Dead" ~ 1,
      vital_status == "Alive" ~ 0,
      TRUE ~ NA_real_  # For "Not Reported" or other cases
    )
  )


z_scored_bones <- as.data.frame(scale(bones))
z_scored_bones$case_submitter_id <- row.names(bones)
combined_data <- merge(merged_data_data_main, z_scored_bones, by = "case_submitter_id", all = FALSE)

head(combined_data)
selected_genes <- colnames(bones)

results <- data.frame(Gene = character(), Hazard_Ratio = numeric(), P_Value = numeric(), stringsAsFactors = FALSE)
survival_data <- combined_data

# Fit Cox model for each gene
for (gene in selected_genes) {
  # Check if the gene exists in survival_data
  if (gene %in% colnames(survival_data)) {
    # Fit Cox model
    cox_model <- coxph(Surv(survival_time, event) ~ survival_data[[gene]], data = survival_data)
    cox_summary <- summary(cox_model)
    
    # Extract hazard ratio and p-value
    hazard_ratio <- cox_summary$coefficients[1, "exp(coef)"]
    p_value <- cox_summary$coefficients[1, "Pr(>|z|)"]
    
    # Store results
    results <- rbind(results, data.frame(Gene = gene, Hazard_Ratio = hazard_ratio, P_Value = p_value))
  }
}
results$Adjusted_P_Value <- p.adjust(results$P_Value, method = "BH")
print(head(results))
significant_genes <- results %>% filter(Adjusted_P_Value < 0.05)
significant_genes$Risk <- ifelse(significant_genes$Hazard_Ratio > 1, "Risky", "Protective")

# Create the bar plot
ggplot(significant_genes, aes(x = reorder(Gene, Hazard_Ratio), y = Hazard_Ratio, fill = Risk)) +
  geom_bar(stat = "identity", width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c("Risky" = "red", "Protective" = "blue")) +
  labs(
    title = "Hazard Ratios of Selected RFE Genes",
    x = "Gene",
    y = "Hazard Ratio",
    fill = "Type"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

srvival_time_index <- which(names(survival_data) == "survival_time")

library(randomForestSRC)

# Subset the data starting from the survival_time column
rf_data <- survival_data[, srvival_time_index:ncol(survival_data)]

rf_data <- na.omit(rf_data)

# Fit random survival forest
rsf_model <- rfsrc(Surv(survival_time, event) ~ ., data = rf_data)

# Print model summary
print(rsf_model)

# Extract variable importance
vimp_values <- vimp(rsf_model)

# Sort variables by importance
sorted_vimp <- sort(vimp_values$importance, decreasing = TRUE)

# Print top variables
print(sorted_vimp)

# Plot top variables
par(mar = c(5, 12, 4, 2))  # Increase left margin (second value)

sorted_vimp_reversed <- rev(sorted_vimp[1:15])  # Reverse the order

# Adjust margins to fit long labels
par(mar = c(5, 12, 4, 2))  # Increase left margin for labels

# Create the horizontal barplot
barplot(
  sorted_vimp_reversed,
  horiz = TRUE,                # Horizontal bars
  las = 1,                     # Keep labels horizontal for readability
  main = "Top 15 Variable Random Forest Importances",
  col = "steelblue",
  cex.names = 0.8              # Adjust label size
)

## genes matching both groups: 
vimp_genes <- names(sorted_vimp_reversed)  # Extract gene names from sorted_vimp_reversed
significant_gene_names <- significant_genes$Gene  # Extract gene names from significant_genes
# Find the intersection of the two groups
overlapping_genes <- intersect(vimp_genes, significant_gene_names)
# Display the overlapping genes
print(overlapping_genes)
