# Load the package
library(minfi)

#####################################################################################################

## location of idat files
idat_files <- list.files("D:/OMIX002042-01.wjz.116s/rawdata.wjz.116s", pattern = "*.idat", full.names = TRUE)
head(idat_files)

## format targets, it is important to just have sample name otherwise minfi thinks there is duplicates due to there be a grn and red per sample
## it looks for green and red per each sample on its own 

targets <- data.frame(
  Basename = gsub("_(Grn|Red)\\.idat", "", idat_files),  # Remove _Grn.idat or _Red.idat
  stringsAsFactors = FALSE
)

targets <- (unique(targets))

## give it unique of column name 
data <- read.metharray.exp(targets = targets, extended = TRUE, verbose = FALSE)

#########################################################################################################
## this requires manifests for the arrays to be installed 

## quality report 
qcReport(data, pdf = "QCReport.pdf")

## filter probes 
detP <- detectionP(data)
length(detP)
keep <- rowMeans(detP <= 0.01) > 0.99
length(keep)

data_filtered <- data[keep, ]

## normalization and filtering 
data_normalized <- preprocessFunnorm(data_filtered)

save(data_normalized, file = "C:/Users/u/Dropbox/SharedDesktopFiles/ShannghaiOsteoMethylation.rds")
