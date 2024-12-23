
library(ggplot2)
library(matrixStats)
library(survival)
library("survminer")
library(data.table)
library(ggplot2)
library(dplyr)
library(reshape)

survival <- read.table("C:/Users/Waleed/Downloads/wcdt_fu_data_231128.txt", sep = " ", header = TRUE)
#samples to keep 

#########################################
#### load needed data 
x <- read.table("C:/Users/Waleed/Downloads/site_wcdt.csv", sep =",", header = TRUE)

xa <- read.table("C:/Users/Waleed/Dropbox/SharedDesktopFiles/Lundberg_at_al/counts_and_TPMs/Counts_TPM_Genes.txt", header =TRUE, check.names = FALSE)
xb <- read.table("C:/Users/Waleed/Dropbox/SharedDesktopFiles/Lundberg_at_al/counts_and_TPMs/Counts_TPM_Genes_PROMOTE.txt", header =TRUE, check.names = FALSE)

x1 <- merge(xa[1:(ncol(xa)-1)],xb[1:(ncol(xb))],by="FEATURE_ID")
rm(xa,xb)

#stemness we chose to look at median z score of SALL4, NANOG, OCT4a (POU5F1) and ALDH1A1 (ALDH1)
### subset data for needed genes
#genes <- c('KLK3','KLK2', 'FKBP5', 'STEAP1', 'STEAP2', 'PLPP1', 'RAB3B', 'NKX3-1', 'ACSL3')
#genes <- c('KLK3','NKX3-1', 'CHRNA2','SLC45A3','TRGC1','TRGC1','NAP1L2', 'KLK2', 'FKBP5', 'STEAP1', 'STEAP2', 'PLPP1', 'RAB3B','ACSL3')
genes <- c('SYP', 'CHGA', 'ENO2')
#genes <- read.table("C:/Users/Waleed/Dropbox/SharedDesktopFiles/NEPCGenesUP_PMID29132337.csv", sep = ",")
#genes <- read.table("C:/Users/Waleed/Dropbox/SharedDesktopFiles/NEPC genes.txt", sep = " ")
#genes <- genes$V1
#genes[89] <- 'IL1B'

x2 <- x1[x1$gene_name %in% genes,]
t1 <- data.frame(t(x2))
colnames(t1) <- x2$gene_name

Zscore = function(x,output){
  m <-  median(as.numeric(x))
  mad <-  mad(as.numeric(x))
  
  if (mad == 0) {
    mad = .0001
  }
  
  n <- vector(mode="numeric", length=length(x))
  for (i in 1:length(x)){
    n[i] <-  (0.6745*(x[i]-m)/ median(mad))
  }
  # return product
  return(n)
}

n.exp <- x2[,2:(ncol(x2)-1)]
n.exp <- log(n.exp+1)
rownames(n.exp) <- x2$gene_name
zexp <- data.frame(apply(n.exp,1,Zscore), check.names = FALSE);
#zexp <- data.frame(t(n.exp)) ## do this if no z scoring
rownames(zexp) <- colnames(n.exp)

## remove NEPC
zexp$sampleID <- paste(rownames(zexp))
zexp <- zexp %>%
  mutate(PatientID = sub("^(\\w+-\\d+)-.*", "\\1", sampleID))

## 
filtered_df <- survival %>% filter(Subtypes == "AR-/NE+" | Subtypes == "AR+/NE+")

## subset 
x_NEPC <- zexp[zexp$PatientID %in% filtered_df$PATIENT_ID,]
x_NEPC <- x_NEPC[,1:(ncol(x_NEPC)-2)]
x_NEPC$NEPC_status <- 'NEPC'

x_NEPCremoved <- zexp[!zexp$PatientID %in% filtered_df$PATIENT_ID,]
x_NEPCremoved <- x_NEPCremoved[,1:(ncol(x_NEPCremoved)-2)]
x_NEPCremoved$NEPC_status <- 'non_NEPC'

## ML
if (!requireNamespace("caret", quietly = TRUE)) install.packages("caret")
if (!requireNamespace("e1071", quietly = TRUE)) install.packages("e1071")
if (!requireNamespace("kernlab", quietly = TRUE)) install.packages("kernlab")
if (!requireNamespace("randomForest", quietly = TRUE)) install.packages("randomForest")

if (!requireNamespace("ROSE", quietly = TRUE)) install.packages("ROSE")

library(caret)
library(e1071)
library(kernlab)
library(randomForest)
library(ROSE)

data <- rbind(x_NEPCremoved, x_NEPC)
#data <- select(data, -c(PatientID, sampleID))

# Convert 'ARstatus' to a factor with descriptive level names
data$NEPC_status <- factor(ifelse(data$NEPC_status == "NEPC", "Positive", "Negative"))

data$NEPC_status <- as.factor(data$NEPC_status)

# Apply SMOTE
smote_result <- ROSE(NEPC_status ~ ., data = data, seed = 1)$data

##
set.seed(123)

# Define cross-validation method
trainControl <- trainControl(method = "cv",
                             number = 10,  # 10-fold cross-validation
                             classProbs = TRUE,  # For ROC/AUC
                             summaryFunction = twoClassSummary)  # For binary classification metrics

# Define a grid of hyperparameters to search over
svmGrid <- expand.grid(sigma = c(0.1, 0.01),  # You may need to adjust these based on your specific dataset
                       C = c(0.5, 1, 2))  # Regularization parameter

# Train the model using cross-validation with SVM
model_cv_svm <- train(NEPC_status ~ ., 
                      data = data, 
                      method = "svmRadial",  # Use SVM with radial basis function kernel
                      tuneGrid = svmGrid,  # Use the defined grid of hyperparameters
                      trControl = trainControl,
                      preProcess = c("center", "scale"),  # Standardize features
                      metric = "ROC")  # Optimize for ROC AUC in this example

# Print the summary of the trained model
print(model_cv_svm)

# Assuming 'testData' is already prepared and 'ARstatus' is a binary numeric variable
# Convert 'ARstatus' in 'testData' to a factor with levels matching the training data
#testData$ARstatus <- factor(testData$ARstatus, levels = c("0", "1"), labels = c("Negative", "Positive"))

# Make predictions using the optimal SVM model

# Splitting data into training and testing sets
set.seed(123)  # For reproducibility
splitIndex <- createDataPartition(data$NEPC_status, p = 0.8, list = FALSE)
trainingData <- data[splitIndex, ]
testingData <- data[-splitIndex, ]


# Predicting on the testing data
predictedClasses <- predict(model_cv_svm, newdata = testingData)

# Creating the confusion matrix
confMat <- confusionMatrix(predictedClasses, testingData$NEPC_status)
print(confMat)

rfModel <- randomForest(NEPC_status ~ ., data=data, importance=TRUE)

# Viewing feature importance
importance(rfModel)

# Plotting feature importance
varImpPlot(rfModel)

