rm(list=ls())

## format and align th emethylation and gene expression of TARGET samples
## combined bone gene expression 
geneEXP <- readRDS("C:/Users/u/Dropbox/SharedDesktopFiles/Bone_CombinedGeneExpression.rds")

methylation <- read.delim("C:/Users/u/Dropbox/TARGET/mergedMethylation.tsv", check.names = FALSE)
print(colnames(methylation))

methylation <- methylation[,1:ncol(methylation)-1] # remove last column it is empty

# View updated column names
print(colnames(beta_matrix))

# update sample names 

samples <- read.delim("C:/Users/u/Dropbox/TARGET/TARGET_methylationSample sheet.tsv", check.names = FALSE)

methylation_colnames <- colnames(methylation)
samples$`File Name` <- gsub("\\.sesame\\.level3betas\\.txt$", "", samples$`File Name`)

# Match column names to 'File Name' and replace with 'Case ID'
updated_colnames <- sapply(methylation_colnames, function(colname) {
  match_index <- match(colname, samples$`File Name`)
  if (!is.na(match_index)) {
    return(samples$`Case ID`[match_index])
  } else {
    return(colname) # Retain original if no match found
  }
})

# Update column names of beta_matrix
colnames(methylation) <- updated_colnames

## select TARGET samples in gene Expression 
geneEXP_sub <- geneEXP[, colnames(geneEXP) %in% colnames(methylation), drop = FALSE]

# Subset methylation to only include columns matching column names in geneEXP
methylation_sub <- methylation[, colnames(methylation) %in% colnames(geneEXP), drop = FALSE]

geneEXP_sub <- cbind(geneEXP[, 1, drop = FALSE], geneEXP_sub)
methylation_sub <- cbind(methylation[, 1, drop = FALSE], methylation_sub)

geneEXP_sub <- geneEXP_sub[, c(1, match(colnames(methylation_sub)[-1], colnames(geneEXP_sub)[-1]) + 1), drop = FALSE]
methylation_sub <- methylation_sub[, c(1, match(colnames(geneEXP_sub)[-1], colnames(methylation_sub)[-1]) + 1), drop = FALSE]

saveRDS(geneEXP_sub, "C:/Users/u/Dropbox/TARGET/geneEXP_matched.rds")

# Save methylation_sub as an .rds file
saveRDS(methylation_sub, "C:/Users/u/Dropbox/TARGET/methylation_matched.rds")

## Extract gene Coordinates and CpGs of genes of interest
library(dplyr)
library(data.table)
library(stringr)
## gene Coordinates 

genes <- c("ITGA10", "PTPRZ1", "IBSP", "MRC2", "SLC8A3", "COL22A1", "COL24A1","COL13A1")

geneCoords <- fread(
  "C:/Users/u/Dropbox/TARGET/gencode.v36.basic.annotation.gtf/gencode.v36.basic.annotation.gtf", 
  sep = "\t", 
  header = FALSE, 
  skip = "#"  # Skips lines starting with '#'
)
names(geneCoords)

geneCoords_filtered <- geneCoords %>%
  filter(V3 == "gene") %>%
  mutate(
    gene_name = str_extract(V9, 'gene_name "([^"]+)"'),
    gene_name = str_replace_all(gene_name, 'gene_name |"', "") # Clean extracted names
  )


# Subset to only include specified genes

genes <- geneCoords_filtered %>%
  filter(gene_name %in% genes) %>%
  select(chromosome = V1, start = V4, end = V5, strand = V7, gene_name)

# View selected genes with genomic coordinates
print(genes)

# Define ±15 kb boundaries
genes <- genes %>%
  mutate(start_15kb = pmax(start - 20000, 0), # Avoid negative start positions
         end_15kb = end + 20000)

# probe manifest 
probes <- read.delim("C:/Users/u/Dropbox/TARGET/HM450.hg38.manifest.gencode.v36.tsv/HM450.hg38.manifest.gencode.v36.tsv")
names(probes)

cpgs_near_genes <- probes %>%
  filter(CpG_chrm %in% genes$chromosome) %>% # Use 'CpG_chrm' for chromosomes
  rowwise() %>%
  filter(any(
    between(CpG_beg, 
            genes$start_15kb[genes$chromosome == CpG_chrm], 
            genes$end_15kb[genes$chromosome == CpG_chrm])
  ))  

## Now subset the gene Expression to only include the slected genes
geneEXP_sub <- geneEXP_sub[geneEXP_sub$gene_id %in% genes$gene_name,]

methylation_sub <- methylation_sub[methylation_sub$CpG_ID %in% cpgs_near_genes$probeID, ]

correlation_results <- data.frame(
  CpG_ID = character(),
  Gene = character(),
  Correlation = numeric(),
  PValue = numeric(),
  stringsAsFactors = FALSE
)

testgene <- genes[8]
Cpg_sub <- cpgs_near_genes[cpgs_near_genes$CpG_chrm == testgene$chromosome, ]
Cpg_sub <- Cpg_sub[Cpg_sub$CpG_beg >= testgene$start,]
Cpg_sub <- Cpg_sub[Cpg_sub$CpG_end <= testgene$end,]

geneEXP_sub_tmp <- geneEXP_sub[geneEXP_sub$gene_id %in% testgene$gene_name ,]
methylation_tmp <- methylation_sub[methylation_sub$CpG_ID %in% Cpg_sub$probeID,  ]

for (j in 1:nrow(methylation_tmp)) {
  cpg_id <- methylation_tmp$CpG_ID[j]
  cpg_values <- as.numeric(methylation_tmp[j, -1]) # Exclude CpG_ID
  gene_expression <- as.numeric(geneEXP_sub_tmp[geneEXP_sub_tmp$gene_id == testgene$gene_name, -1]) # Exclude gene_id
  
  # Remove NA values by matching non-NA indices
  complete_cases <- complete.cases(cpg_values, gene_expression)
  cpg_values <- cpg_values[complete_cases]
  gene_expression <- gene_expression[complete_cases]
  
  # Check if any values remain after filtering
  if (length(cpg_values) > 0 && length(gene_expression) > 0) {
    # Perform correlation test
    correlation_test <- cor.test(cpg_values, gene_expression, method = "pearson")
    correlation_results <- rbind(correlation_results, data.frame(
      CpG_ID = cpg_id,
      Gene = testgene$gene_name,
      Correlation = correlation_test$estimate,
      PValue = correlation_test$p.value
    ))
  } else {
    cat("No matching non-NA values for CpG:", cpg_id, "and gene:", testgene$gene_name, "\n")
  }
}

print(correlation_results)

write.table(correlation_results, "C:/Users/u/Dropbox/SharedDesktopFiles/MethylationCorr.Results", col.names = TRUE, row.names = FALSE)

sig <- correlation_results[correlation_results$PValue < 0.05,]

sub <- sig[sig$Correlation < 0,  ] ## change this to select positive or negative correlations

sub$Correlation <- as.numeric(sub$Correlation)
sub <- sub %>%
  arrange(sub$Correlation)

cat(paste(round(sub$Correlation, digits =2), collapse = "\n"))
cat(paste(sub$CpG_ID, collapse = "\n"))

#cat(paste(round(sub$adj.P, digits =15), collapse = "\n"))
cat(paste(format(sub$PValue, scientific = TRUE, digits = 3), collapse = "\n"))


write.table(correlation_results, "C:/Users/u/Dropbox/SharedDesktopFiles/MethylationCorr.SigResults", col.names = TRUE, row.names = FALSE)

# View the correlation results
print(correlation_results)

## Annotate CpGS

library(GenomicRanges)
library(AnnotationHub)
library(rtracklayer)
library(TFBSTools)
library(Biostrings)
library(GenomicRanges)
library(JASPAR2022)
library(BSgenome.Hsapiens.UCSC.hg38)
genome <- BSgenome.Hsapiens.UCSC.hg38

sub <- sig[sig$Correlation > 0,  ] ## change this to select positive or negative correlations

## select original CpG object with CpG coordinate info 
# subset it to only include CpGs in significant subset (positive or negative) 

CpG_subset <- cpgs_near_genes[cpgs_near_genes$probeID %in% sub$CpG_ID,]
head(CpG_subset)

# make a GRanges Object, needed for Downstream analysis

cpg_gr <- GRanges(
  seqnames = CpG_subset$CpG_chrm,
  ranges = IRanges(start = CpG_subset$CpG_beg, end = CpG_subset$CpG_end),
  strand = CpG_subset$probe_strand,
  probeID = CpG_subset$probeID,
  genesUniq = CpG_subset$genesUniq,
  geneNames = CpG_subset$geneNames,
  transcriptTypes = CpG_subset$transcriptTypes,
  transcriptIDs = CpG_subset$transcriptIDs,
  distToTSS = CpG_subset$distToTSS
)

# Encode cis regulatory element file downloaded doi:10.17989/ENCSR439EAZ
## https://www.encodeproject.org/annotations/ENCSR439EAZ/

# JASPAR vertebrate TF PWM dowloaded
# https://jaspar2022.genereg.net/downloads/

### Enhancer annotation
raw_data <- fread("F:/Annotations/ENCODE_Enhancers/ENCFF924IMH.bed", header = FALSE)

# Keep only the first three columns (chr, start, end)
bed_data <- raw_data[, .(chr = V1, start = V2, end = V3)]

# Convert to GRanges
encode_cre <- GRanges(
  seqnames = bed_data$chr,
  ranges = IRanges(start = bed_data$start, end = bed_data$end)
)

# 10th col contains info of type of Enhancer 
encode_cre$annotation <- raw_data$V10

## expand CPG region and get corresponding sequence
# CPG region expanded 1000 (+/-500) to look for nearby enhancers

expanded_cpg_gr <- resize(cpg_gr, width = width(cpg_gr) + 1000, fix = "center")

# GRAnges function to ensure it is within chromosome boundaries
expanded_cpg_gr <- trim(expanded_cpg_gr)

expanded_cpg_sequences <- getSeq(genome, expanded_cpg_gr)
overlaps <- findOverlaps(cpg_gr, expanded_cpg_gr)

# Annotate CpGs with overlap information
cpg_gr$Enhancer <- ifelse(seq_along(cpg_gr) %in% queryHits(overlaps), "Yes", "No")

# Initialize annotation column for CpGs
cpg_gr$Enhancer_Annotation <- NA

# Map enhancer annotations to overlapping CpGs
overlap_indices <- queryHits(overlaps)
cpg_gr$Enhancer_Annotation[overlap_indices] <- encode_cre$annotation[subjectHits(overlaps)]

annotated_cpg <- as.data.frame(cpg_gr)

## 

pfm_list <- getMatrixSet(
  JASPAR2022,
  opts = list(collection = "CORE", tax_group = "vertebrates", matrixtype = "PFM")
)

tf_names <- c("RUNX2", "SOX9", "SP7", "TWIST1", "ETS1", "NFATC1", 
              "MYC", "FOXO3", "MEF2C", "SMAD1", "GLI2", "HOXA10", 
              "ZEB1", "CREB1", "TCF7")

# Retrieve PFMs from JASPAR
pfm_list <- lapply(tf_names, function(tf_name) {
  query <- list(name = tf_name, tax_group = "vertebrates", collection = "CORE")
  tryCatch(getMatrixSet(JASPAR2022, opts = query)[[1]], error = function(e) NULL)
})

# Filter out NULLs (in case a TF PFM is not found)
names(pfm_list) <- tf_names
pfm_list <- Filter(Negate(is.null), pfm_list)

# Convert PFMs to PWMs
pwm_list <- lapply(pfm_list, toPWM)


tf_binding_results <- list()

# Scan for each TF
for (tf_name in names(pwm_list)) {
  pwm <- pwm_list[[tf_name]]
  
  # Match PWM to expanded CpG sequences
  matches <- lapply(expanded_cpg_sequences, function(seq) {
    matchPWM(pwm@profileMatrix, seq, min.score = "70%")  # Adjust min.score as needed
  })
  
  # Annotate CpGs with the number of matches
  tf_binding_results[[tf_name]] <- sapply(matches, length)
}
total_binding_sites <- sapply(tf_binding_results, sum)
print(total_binding_sites)

binding_proportions <- sapply(tf_binding_results, function(x) mean(x > 0))
print(binding_proportions)

multi_tf_cpgs <- rowSums(sapply(tf_binding_results, function(x) x > 0)) > 1
print(which(multi_tf_cpgs))  # Indices of CpGs with multiple TF matches

cpg_gr$TFs_Binding <- apply(sapply(tf_binding_results, function(x) ifelse(x > 0, names(tf_binding_results), "")), 1, function(x) paste(unique(x[x != ""]), collapse = ";"))

annotated_cpg <- as.data.frame(cpg_gr)
names(annotated_cpg)

## Subset only CPGS that contain TFB or enhancer annotation
cpg_with_enhancer_or_tfbs <- subset(
  annotated_cpg,
  Enhancer == "Yes" | TFs_Binding != ""
)

names(sub)
Annotated_sub <- sub[sub$CpG_ID %in% cpg_with_enhancer_or_tfbs$probeID ,]

Annotated_sub <- merge(
  Annotated_sub,
  cpg_with_enhancer_or_tfbs[, c("probeID", "Enhancer_Annotation", "TFs_Binding")],
  by.x = "CpG_ID",
  by.y = "probeID",
  all.x = TRUE
)

# Create a count summary for Enhancer and TFBS
Annotated_sub$Enhancer_Status <- ifelse(Annotated_sub$Enhancer_Annotation != "", "Enhancer", "No Enhancer")
Annotated_sub$TFBS_Status <- ifelse(Annotated_sub$TFs_Binding != "", "TFBS", "No TFBS")

Annotated_sub$Correlation <- as.numeric(Annotated_sub$Correlation)

Annotated_sub <- Annotated_sub %>%
  arrange(desc(Annotated_sub$Correlation))

cat(paste(round(Annotated_sub$Correlation, digits =2), collapse = "\n"))

cat(paste(Annotated_sub$Enhancer_Annotation, collapse = "\n"))

#cat(paste(round(sub$adj.P, digits =15), collapse = "\n"))
cat(paste(format(Annotated_sub$PValue, scientific = TRUE, digits = 3), collapse = "\n"))

# View the resulting subset
head(cpg_with_enhancer_or_tfbs)
