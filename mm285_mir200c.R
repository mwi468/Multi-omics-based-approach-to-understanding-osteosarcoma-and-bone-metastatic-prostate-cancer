rm(list = ls())

# Load necessary libraries
library(sesame)
library(ggplot2)
library(wheatmap)
library(readxl)

# Load beta values and color data
betass <- readRDS("/Users/iqbalw/Downloads/20210319_1656_mouse_array_samples_NAchanged.rds")
color_df <- readRDS("/Users/iqbalw/Downloads/20210324_mouse_array_samples_GroupColorCodes.rds")

# Subset betass by available sample IDs
mat_keep_rows <- as.list(color_df$Sample_ID)
betass <- betass[, colnames(betass) %in% mat_keep_rows]

# Update color_df to match the remaining samples in betass
mat_keep_rows2 <- as.list(colnames(betass))
color_df <- color_df[color_df$Sample_ID %in% mat_keep_rows2, ]

# Order by Sample ID in both betass and color_df
betass <- betass[, order(colnames(betass))]
color_df <- color_df[order(color_df$Sample_ID), ]

# Initialize an empty data frame for distance and probe count
df <- data.frame(matrix(vector(), 0, 2, dimnames = list(c(), c("Distance", "Probe_Count"))), stringsAsFactors = FALSE)

# Define distance range and loop through to visualize gene data
dist <- 1500
while (dist <= 20000) {
  b1 <- visualizeGene('Mir200c', betass, refversion = "mm10", platform = "MM285", 
                      upstream = dist, dwstream = dist, cluster.samples = TRUE, draw = FALSE)
  x <- data.frame(c(dist, dim(b1)[1]))
  x <- t(x)
  colnames(x) <- c("Distance", "Probe_Count")
  df <- rbind(df, x)
  dist <- dist + 200
}

# Initialize a second data frame to get distances with changes in probe count
df2 <- data.frame(matrix(vector(), 0, 2, dimnames = list(c(), c("Distance", "Probe_Count"))), stringsAsFactors = FALSE)

# Get unique probe counts and their max distances
j <- unique(df$Probe_Count)
for (i in 1:length(j)) {
  d <- df[df$Probe_Count == j[i], ]
  di <- max(d$Distance)
  xi <- data.frame(c(j[i], di))
  xi <- t(xi)
  colnames(xi) <- c("Probe_Count", "Distance")
  df2 <- rbind(df2, xi)
}

# Load additional metadata
tissue_color <- read_excel('~/samplesheets/2021/20210222_Color_Scheme.xlsx') %>% 
  with(setNames(HEX, TISSUE))
meta <- read_excel("~/samplesheets/2021/20210226_MouseArray_SampleTableV4.xlsx")

# Prepare for distance-based visualization
dist <- unique(df2$Distance)
test <- data.frame(matrix(vector(), ncol(betass), 0), stringsAsFactors = FALSE)

# Main distance for final visualization
distance <- 20100
b1 <- visualizeGene('Mir200c', betass, refversion = "mm10", platform = "MM285", 
                    upstream = distance, dwstream = distance, cluster.samples = FALSE, draw = FALSE)
b1 <- data.frame(t(b1))

# Calculate median for each probe
b1$ProbeMedian <- apply(b1, 1, median)

# Add tissue information
b1$Sample_ID <- rownames(b1)
sub <- color_df[, 1:2]
b1 <- merge(b1, sub, by.y = 'Sample_ID', sort = FALSE)

# Loop through unique tissues and order by Probe Median
x <- unique(b1$TISSUE)
s <- NULL

for (i in 1:length(x)) {
  if (i == 1) {
    s <- b1[b1$TISSUE == x[i], ]
    s <- s[order(s$ProbeMedian, decreasing = TRUE), ]
    s$TISSUE_Median <- median(as.numeric(s$ProbeMedian), na.rm = TRUE)
  } else {
    s2 <- b1[b1$TISSUE == x[i], ]
    s2 <- s2[order(s2$ProbeMedian, decreasing = TRUE), ]
    s2$TISSUE_Median <- median(as.numeric(s2$ProbeMedian), na.rm = TRUE)
    s <- rbind(s, s2)
  }
}

# Order tissues based on their median probe value
s <- s[order(s$TISSUE_Median, decreasing = TRUE), ]
s$MiR200_MeanMeth <- apply(s[, (grep("cg43961439_BC21", names(s)):grep("cg43961464_TC11", names(s)))], 1, mean)
test <- cbind(test, s$MiR200_MeanMeth)

# Filter out NA tissues
s <- s[s$TISSUE != 'NA', ]

# Rearrange so rows are samples
rownames(s) <- s$Sample_ID
fs <- s[, 2:(ncol(s) - 4)]
fs <- data.frame(t(fs), check.names = FALSE)

# Generate the heatmap
x <- visualizeGene('Mir200c', fs, refversion = "mm10", platform = "MM285", 
                   show.sampleNames = FALSE, upstream = distance, dwstream = distance, 
                   cluster.samples = FALSE) +
  WColorBarV(meta$Tissue_Corrected[match(colnames(fs), meta$IDAT)], 
              cmp = CMPar(label2color = tissue_color), RightOf('betas')) +
  WLegendV(NULL, BottomRightOf('betas', h.pad = 0.2), height = 0.06) +
  WColorBarV(s$MiR200_MeanMeth, LeftOf('betas')) + 
  WLegendV(NULL, name = "", BottomLeftOf('betas', h.pad = -0.4), height = 0.06)

# Save the heatmap as a PNG
png(paste("~/gallery/20210614_MIR_ManualSorting", distance, ".png"), height = 20, width = 10, units = 'in', res = 300)
print(x)
dev.off()

# List of probes specific to Mir200c (for reference)
# "cg43961439_BC21", "cg43961440_BC21", "cg43961441_BC21", "cg43961451_BC21", "cg43961464_TC11"

