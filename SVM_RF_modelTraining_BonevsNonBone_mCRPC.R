rm(list=ls())
library(tidyr)
library(dplyr)
library(reshape2)
library(RColorBrewer)

######## Model Training

#### load combined data and annotation files 
Combined_Annotations  <- read.delim("C:/Users/u/Dropbox/SharedDesktopFiles/Datasets/ProstateCancer_Human/CombinedAnnotations.tsv", header =TRUE, check.names = FALSE)
numeric_data_matrix <- read.delim("C:/Users/u/Dropbox/SharedDesktopFiles/Datasets/ProstateCancer_Human/CombinedNumerical_data.tsv", header =TRUE, check.names = FALSE)

## Filter tissue types
unwanted_sites <- c("Other", "unknown", "Primary", "Prostate")

filtered_data <- Combined_Annotations %>%
  filter(!TISSUE_SITE %in% unwanted_sites)

## filter the numerical matrix based on this: 
sample_ids_to_keep <- filtered_data$SAMPLE_ID

# Filter columns in `numeric_data_matrix` to include only the sample IDs in `sample_ids_to_keep`
filtered_matrix <- numeric_data_matrix[rownames(numeric_data_matrix) %in% sample_ids_to_keep,]

## use SVM, RFE )recursive feature elimnation  and Random Forest to differentiate between bone and non bone 


################## svm 

if (!require("e1071")) install.packages("e1071")
library(e1071)
library(caret)

#############

filtered_data$LABEL <- ifelse(filtered_data$TISSUE_SITE == "Bone", 1, 0)

# Ensure filtered_data only includes samples present in filtered_matrix
filtered_data <- filtered_data[filtered_data$SAMPLE_ID %in% rownames(filtered_matrix), ]
filtered_matrix <- filtered_matrix[rownames(filtered_matrix) %in% filtered_data$SAMPLE_ID, ]

# Reorder filtered_data to match the sample order in filtered_matrix
filtered_data <- filtered_data[match(rownames(filtered_matrix), filtered_data$SAMPLE_ID), ]

# Create the label vector 
labels <- as.factor(filtered_data$LABEL)

##### svm
set.seed(42)

train_control <- trainControl(method = "cv", number = 10)  # 10-fold cross-validation

# Train the SVM model using caret with cross-validation
set.seed(42)

## Croiss validation using linear and non linear kernal to select best model

#svm_linear_grid <- expand.grid(C = seq(0.001, 100, by = 1))
#svm_linear_grid <- expand.grid(C = seq(0.0001, 0.01, by = 0.0005))
#svm_rbf_grid <- expand.grid(
#  C = seq(0.1, 10, by = 0.5),
#  sigma = seq(0.01, 1, by = 0.1)
#)

# Train the SVM model using caret with cross-validation
#svm_linear_tuned <- train(
#  x = filtered_matrix,
#  y = labels,
#  method = "svmLinear",
#  trControl = train_control,
#  tuneGrid = svm_linear_grid
#)

#set.seed(42)  # For reproducibility
#svm_rbf_tuned <- train(
#  x = filtered_matrix,
#  y = labels,
#  method = "svmRadial",
#  trControl = train_control,
#  tuneGrid = svm_rbf_grid
#)
# Print the tuned model results
#print(svm_linear_tuned)
#print(svm_rbf_tuned)

## the best model was linaer with tune grid = .0011

set.seed(42)
train_control <- trainControl(
  method = "cv",         # Cross-validation
  number = 10,           # 10-fold cross-validation
  verboseIter = TRUE     # Display progress
)

# Run the linear SVM model with the best parameter (C = 0.0011)

zero_threshold <- 0.75  # Genes with more than 90% zeros will be removed
variance_threshold <- 0.01  # Genes with variance below this value will be removed
# Remove genes with too many zeros
gene_zero_fraction <- rowMeans(filtered_matrix == 0)
filtered_matrix <- filtered_matrix[gene_zero_fraction <= zero_threshold, ]

# Remove genes with very low variance
gene_variances <- apply(filtered_matrix, 1, var)
filtered_matrix <- filtered_matrix[gene_variances > variance_threshold, ]

set.seed(42)

stratified_control <- trainControl(
  method = "repeatedcv",  # Use repeated cross-validation
  number = 10,           # Number of folds
  repeats = 3,           # Number of repeats
  classProbs = TRUE,     # Enable class probabilities
  summaryFunction = twoClassSummary,  # Use metrics like accuracy, sensitivity, and specificity
  savePredictions = "final"           # Save predictions for further analysis
)

# Train SVM with linear kernel using stratified k-fold CV
set.seed(123)  # Set a seed for reproducibility

labels <- factor(labels, levels = c("0", "1"), labels = c("Class0", "Class1"))
svm_stratified <- train(
  x = filtered_matrix,        # Feature matrix
  y = labels,                 # Class labels
  method = "svmLinear",       # SVM with linear kernel
  trControl = stratified_control,  # Stratified cross-validation
  tuneGrid = data.frame(C = 0.00025),  # Specify the best 'C' value from previous tuning
  metric = "ROC"              # Use ROC for model evaluation
)

print(svm_stratified)
# Access resampling results
resampling_results <- svm_stratified$resample
print(resampling_results)


## Random Forest 

set.seed(42)
rf_stratified_model <- train(
  x = filtered_matrix,              # Feature matrix
  y = labels,                       # Labels (bone vs non-bone)
  method = "rf",                    # Random Forest
  metric = "ROC",                   # Optimize for ROC
  trControl = stratified_control    # Use the stratified control object
)

# Print the trained model details
print(rf_stratified_model)

### Save RF and SVM models 
saveRDS(rf_stratified_model, file = "rf_stratified_model.rds")
# Save the SVM model
saveRDS(svm_stratified, file = "svm_stratified_model.rds")
