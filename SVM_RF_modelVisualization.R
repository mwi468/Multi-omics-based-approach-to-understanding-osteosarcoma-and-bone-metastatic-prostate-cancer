rm(list=ls())

rf_stratified <- readRDS("rf_stratified_model.rds")
svm_stratified <- readRDS("svm_stratified_model.rds")


library(pROC)
library(caret)
library(ggplot2)
library(reshape2)

# Filter predictions for SVM
svm_preds <- svm_stratified$pred[complete.cases(svm_stratified$pred), ]

# Filter predictions for RF
rf_preds <- rf_stratified$pred[complete.cases(rf_stratified$pred), ]

######## Generate ROC curves
roc_svm <- roc(
  response = svm_preds$obs,
  predictor = svm_preds$Class1,
  levels = rev(levels(svm_preds$obs))
)
roc_rf <- roc(
  response = rf_preds$obs,
  predictor = rf_preds$Class1,
  levels = rev(levels(rf_preds$obs))
)

# Combine ROC curves on the same plot
ggplot() +
  geom_line(aes(x = roc_svm$specificities, y = roc_svm$sensitivities, color = "SVM")) +
  geom_line(aes(x = roc_rf$specificities, y = roc_rf$sensitivities, color = "Random Forest")) +
  labs(
    title = "ROC Curve Comparison: SVM vs Random Forest",
    x = "1 - Specificity",
    y = "Sensitivity",
    color = "Model"
  ) +
  annotate("text", x = 0.6, y = 0.1, label = paste("SVM AUC =", round(auc(roc_svm), 3)), color = "blue") +
  annotate("text", x = 0.6, y = 0.2, label = paste("RF AUC =", round(auc(roc_rf), 3)), color = "green") +
  theme_minimal()

######### gENERATE cONFUSION mATRIX 
# Add a Model column to distinguish between SVM and Random Forest
svm_cm <- table(Prediction = svm_preds$pred, Reference = svm_preds$obs)
rf_cm <- table(Prediction = rf_preds$pred, Reference = rf_preds$obs)

svm_cm_df <- as.data.frame(svm_cm)
rf_cm_df <- as.data.frame(rf_cm)

svm_cm_df$Model <- "SVM"
rf_cm_df$Model <- "Random Forest"

# Ensure both data frames have the same structure
colnames(svm_cm_df) <- c("Prediction", "Reference", "Frequency", "Model")
colnames(rf_cm_df) <- c("Prediction", "Reference", "Frequency", "Model")

# Combine data for faceting
combined_cm_df <- rbind(svm_cm_df, rf_cm_df)

# Rename the factor levels for better readability
combined_cm_df$Prediction <- factor(combined_cm_df$Prediction, levels = c("Class1", "Class0"), labels = c("Bone", "Non-Bone"))
combined_cm_df$Reference <- factor(combined_cm_df$Reference, levels = c("Class1", "Class0"), labels = c("Bone", "Non-Bone"))

# Plot confusion matrix heatmaps in the same row (two columns)
ggplot(data = combined_cm_df, aes(x = Reference, y = Prediction, fill = Frequency)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Frequency), color = "white", size = 4) +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  facet_grid(~ Model) +  # Arrange plots in one row
  labs(title = "Confusion Matrix Heatmaps: SVM vs. Random Forest",
       x = "True Labels",
       y = "Predicted Labels") +
  theme_minimal() +
  theme(strip.text = element_text(size = 14, face = "bold"))

######## printing overall metrics 

library(ggplot2)

svm_metrics <- svm_stratified$resample %>%
  select(ROC, Sens, Spec) %>%
  mutate(Model = "SVM")  # Add a column for the model name

rf_metrics <- rf_stratified$resample %>%
  select(ROC, Sens, Spec) %>%
  mutate(Model = "Random Forest")  # Add a column for the model name

# Combine the metrics from both models
combined_metrics <- bind_rows(
  svm_metrics %>%
    pivot_longer(cols = c(ROC, Sens, Spec), names_to = "Metric", values_to = "Value"),
  rf_metrics %>%
    pivot_longer(cols = c(ROC, Sens, Spec), names_to = "Metric", values_to = "Value")
)

ggplot(combined_metrics, aes(x = Metric, y = Value, fill = Model)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Set2") +  # Use a nice color palette
  labs(
    title = "Comparison of Metrics: SVM vs. Random Forest",
    x = "Metric",
    y = "Value"
  ) +
  theme_minimal() +
  theme(
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10)
  )
