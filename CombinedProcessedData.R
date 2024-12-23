rm(list=ls())

library(data.table)
library(dplyr)


##### functions to utilize 

######## Function to aggregate transcript TPM to gene level TPM
## assumes first col is gene names, subsequent cols is data for sample (rows are genes, cols are samples)

aggregate_transcript_to_gene_level <- function(data) {
  # Ensure there is at least one gene column and one sample column
  if (ncol(data) < 2) {
    stop("Data must have at least one gene column and one sample column with numeric values.")
  }
  
  # Rename the first column to "gene_id" for consistency
  colnames(data)[1] <- "gene_id"
  
  # Aggregate values by summing across entries for each gene and sample
  gene_level_data <- data %>%
    group_by(gene_id) %>%
    summarise(across(everything(), \(x) sum(x, na.rm = TRUE))) %>%
    ungroup()
  
  return(gene_level_data)

}

######### function to convert other datasets to TPM from FPKM

## assumes first col is gene names, subsequent cols is data for sample (rows are genes, cols are samples)

convert_fpkm_to_tpm <- function(fpkm_data) {
  # Copy the data frame to store TPM values
  tpm_data <- fpkm_data
  
  # Loop through each column, starting from the second column
  for (i in 2:ncol(fpkm_data)) {
    # Sum of FPKM values in the column (sample)
    total_fpkm <- sum(as.numeric(fpkm_data[, i]), na.rm = TRUE)
    
    # Convert FPKM to TPM
    tpm_data[, i] <- (as.numeric(fpkm_data[, i]) / total_fpkm) * 1e6
  }
  
  return(tpm_data)
}

## Function to clean names

clean_names <- function(x) {
  x <- gsub("[^[:alnum:]_]", "", x)   # Remove special characters except underscores
  x <- trimws(x)                      # Remove leading and trailing whitespace
  toupper(x)                           # Convert to uppercase (optional, if needed)
}

############ combine datasets

########################################### Combine Gene Expression Files
#USCF 
xa <- read.table("C:/Users/Waleed/Dropbox/SharedDesktopFiles/Lundberg_at_al/counts_and_TPMs/Counts_TPM_Genes.txt", header =TRUE, check.names = FALSE)
xb <- read.table("C:/Users/Waleed/Dropbox/SharedDesktopFiles/Lundberg_at_al/counts_and_TPMs/Counts_TPM_Genes_PROMOTE.txt", header =TRUE, check.names = FALSE)

USCF <- merge(xa[1:(ncol(xa)-1)], xb[1:(ncol(xb))], by="FEATURE_ID")

## This data has transcript gene expressions, we should keep only canonical transcripts per gene 

# replace ensemble feature id with gene names (last col)
USCF$FEATURE_ID <- USCF[, ncol(USCF)]

# Remove the last column (original gene names column) from the data frame
USCF <- USCF[, -ncol(USCF)]

# Rename the first column to "Gene_Name" if desired
colnames(USCF)[1] <- "Gene_Name"
rm(xa, xb)  # Remove temporary data frames from memory

USCF_TPM <- aggregate_transcript_to_gene_level(USCF)

## SU2C
x <- read.table("C:/Users/Waleed/Dropbox/SharedDesktopFiles/Datasets/ProstateCancer_Human/SU2C/data_mrna_seq_fpkm_capture.txt", header = TRUE)
SUC_TPM <- convert_fpkm_to_tpm(x)
colnames(SUC_TPM)[1] <- "Gene_Name"
SUC_TPM <- aggregate_transcript_to_gene_level(SUC_TPM)

## WCM_NEPC
x <- read.delim("C:/Users/Waleed/Dropbox/SharedDesktopFiles/Datasets/ProstateCancer_Human/WCM_NEPC/data_mrna_seq_fpkm.txt", header = TRUE)
x <- x[, -2] ## keep only Gene names (1st col) remove entrez id (second column)

WCM_TPM <- convert_fpkm_to_tpm(x)
colnames(WCM_TPM)[1] <- "Gene_Name"
WCM_TPM <- aggregate_transcript_to_gene_level(WCM_TPM)

### cobine datasets
combined_data <- merge(USCF_TPM, SUC_TPM, by = "gene_id", all = TRUE)
combined_data <- merge(combined_data, WCM_TPM, by = "gene_id", all = TRUE)

## reomove any samples that might be duplicated (in case diff studies possess the same sample)

duplicate_cols <- duplicated(as.list(combined_data))

# View which columns are duplicates
duplicated_sample_names <- colnames(combined_data)[duplicate_cols]
print(duplicated_sample_names)

combined_data_unique <- test[, !duplicate_cols]

t_combined_data_unique <- data.frame(t(combined_data_unique[,2:ncol(combined_data_unique)]))

colnames(t_combined_data_unique) <- combined_data_unique$gene_id
rownames(t_combined_data_unique) <- clean_names(colnames(combined_data_unique)[2:ncol(combined_data_unique)])

numeric_data <- t_combined_data_unique[, sapply(t_combined_data_unique, is.numeric)]

# Define the file path where you want to save the data
file_path <- "C:/Users/Waleed/Dropbox/SharedDesktopFiles/Datasets/ProstateCancer_Human/CombinedNumerical_data.tsv"
write.table(numeric_data, file = file_path, sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)


############################################ prepare sample annotation file

########## UCSF expression 

# Tissue site info 
UCSFmeta <- read.csv("C:/Users/Waleed/Dropbox/SharedDesktopFiles/Lundberg_at_al/site_wcdt.csv")

# Clean sample names in UCSF meta
UCSFmeta$SAMPLE_ID <- clean_names(UCSFmeta$SAMPLE_ID)

# Create a data frame from the cleaned column names of UCSF_TPM, excluding the first column
sample_info_UCSF <- data.frame(SAMPLE_ID = clean_names(colnames(USCF_TPM)[2:ncol(USCF_TPM)]))

# Merge sample_info_UCSF with UCSFmeta based on SAMPLE_ID to get Site_detailed
sample_info_UCSF <- merge(sample_info_UCSF, UCSFmeta[, c("SAMPLE_ID", "Site_detailed")], 
                          by = "SAMPLE_ID", all.x = TRUE)

## Subtype
UCSFmeta2 <- read.table("C:/Users/Waleed/Dropbox/SharedDesktopFiles/Lundberg_at_al/wcdt_fu_data_231128.txt", header =TRUE)

# Clean sample names in UCSF meta
UCSFmeta2$SAMPLE_ID <- clean_names(UCSFmeta2$SAMPLE_ID)

# Merge sample_info_UCSF with UCSFmeta based on SAMPLE_ID to get Site_detailed
sample_info_UCSF <- merge(sample_info_UCSF, UCSFmeta2[, c("SAMPLE_ID", "Subtypes")], 
                          by = "SAMPLE_ID", all.x = TRUE)

sample_info_UCSF$Subtypes <- ifelse(
  is.na(sample_info_UCSF$Subtypes), 
  "unknown",  # If the value is NA, set to "unknown"
  ifelse(
    sample_info_UCSF$Subtypes %in% c("AR-/NE+", "AR+/NE+"), 
    "NEPC",    # If AR-/NE+ or AR+/NE+, set to "NEPC"
    "NON_NEPC"  # Otherwise, set to "NON_NEPC"
  )
)

colnames(sample_info_UCSF)[colnames(sample_info_UCSF) == "Site_detailed"] <- "TISSUE_SITE"

# Fill 'Tissue_SITE' with "unknown" where there are NA values
if(any(is.na(sample_info_UCSF$TISSUE_SITE))) {
  sample_info_UCSF$TISSUE_SITE[is.na(sample_info_UCSF$TISSUE_SITE)] <- "unknown"
} else {
  message("No NA values in TISSUE_SITE column.")
}

# Add 'Study' and 'Tumor' columns
sample_info_UCSF$Study <- 'UCSF'
sample_info_UCSF$Tumor <- 'ProstateCancerMetastases'

colnames(sample_info_UCSF)[colnames(sample_info_UCSF) == "Site_detailed"] <- "TISSUE_SITE"

################################ SUC
SUC_meta <- read.delim("C:/Users/Waleed/Dropbox/SharedDesktopFiles/Datasets/ProstateCancer_Human/SU2C/data_clinical_sample.txt", 
                       header = TRUE, skip = 4)

SUC_meta$SAMPLE_ID <- clean_names(SUC_meta$SAMPLE_ID)

# Create a data frame from the cleaned column names of SUC_TPM, excluding the first column
sample_info_SUC <- data.frame(SAMPLE_ID = clean_names(colnames(SUC_TPM)[2:ncol(SUC_TPM)]))

# Merge sample_info_SUC with SUC_meta based on SAMPLE_ID to get the TISSUE_SITE column
sample_info_SUC <- merge(sample_info_SUC, SUC_meta[, c("SAMPLE_ID", "TISSUE_SITE", "NEUROENDOCRINE_FEATURES")], 
                         by = "SAMPLE_ID", all.x = TRUE)

colnames(sample_info_SUC)[colnames(sample_info_SUC) == "NEUROENDOCRINE_FEATURES"] <- "Subtypes"

sample_info_SUC$Subtypes <- ifelse(
  is.na(sample_info_SUC$Subtypes), 
  "unknown",  # If the value is NA, set to "unknown"
  ifelse(
    sample_info_SUC$Subtypes %in% c("Yes"), 
    "NEPC",    # If Yes, set to "NEPC"
    "NON_NEPC"  # Otherwise, set to "NON_NEPC"
  )
)

# Fill in missing TISSUE_SITE values with "unknown" where NA
if(any(is.na(sample_info_SUC$TISSUE_SITE))) {
  sample_info_SUC$TISSUE_SITE[is.na(sample_info_SUC$TISSUE_SITE )] <- "unknown"
} else {
  message("No NA values in TISSUE_SITE column.")
}

## Add study and tumor Info
sample_info_SUC$Study <- 'SUC'
sample_info_SUC$Tumor <- 'ProstateCancerMetastases'

############################################# WCM annotation 
WCM_meta <- read.delim("C:/Users/Waleed/Dropbox/SharedDesktopFiles/Datasets/ProstateCancer_Human/WCM_NEPC/data_clinical_sample.txt", 
                                   header = TRUE, skip = 4)

WCM_meta$SAMPLE_ID <- clean_names(WCM_meta$SAMPLE_ID)

sample_info_WCM <- data.frame(SAMPLE_ID = clean_names(colnames(WCM_TPM)[2:ncol(WCM_TPM)]))

# Merge sample_info_WCM with WCM_meta based on SAMPLE_ID to get Tissue Site information
# Assuming the tissue site column in WCM_meta is called "TISSUE_SITE"
sample_info_WCM <- merge(sample_info_WCM, WCM_meta[, c("SAMPLE_ID", "TUMOR_TISSUE_SITE","DISEASE_CODE" )], 
                         by = "SAMPLE_ID", all.x = TRUE)

colnames(sample_info_WCM)[colnames(sample_info_WCM) == "TUMOR_TISSUE_SITE"] <- "TISSUE_SITE"
colnames(sample_info_WCM)[colnames(sample_info_WCM) == "DISEASE_CODE"] <- "Subtypes"

### clean up NEPC info 

sample_info_WCM$Subtypes <- ifelse(
  is.na(sample_info_WCM$Subtypes), 
  "unknown",  # If the value is NA, set to "unknown"
  ifelse(
    sample_info_WCM$Subtypes %in% c("CRPC-NE"), 
    "NEPC",    # CRPC-NE 
    "NON_NEPC"  # Otherwise, set to "NON_NEPC"
  )
)



# Fill 'TISSUE_SITE' with "unknown" where there are NA values
sample_info_WCM$TISSUE_SITE[is.na(sample_info_WCM$TISSUE_SITE)] <- "unknown"

# Add 'Study' and 'Tumor' columns
sample_info_WCM$Study <- 'WCM'
sample_info_WCM$Tumor <- 'ProstateCancer'

#########################################################################################
## Combine meta data from all studies

Combined_Annotations <- rbind(sample_info_WCM,sample_info_UCSF,sample_info_SUC)
Combined_Annotations$TISSUE_SITE <- gsub("^LN$|^Lymph node$", "Lymph Node", Combined_Annotations$TISSUE_SITE, ignore.case = TRUE)



#### save Annotation file 
file_path <- "C:/Users/Waleed/Dropbox/SharedDesktopFiles/Datasets/ProstateCancer_Human/CombinedAnnotations.tsv"
write.table(Combined_Annotations, file = file_path, sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
