# Clear the environment
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

# Visualize gene data for Igf2
b1 <- visualizeGene('Igf2', betass, refversion = "mm10", platform = "MM285", 
                    cluster.samples = TRUE, dwstream = 80000, draw = FALSE)

# Transpose the data frame and calculate median for each probe
b1 <- data.frame(t(b1))
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

# Order the tissues based on their median probe value
s <- s[order(s$TISSUE_Median, decreasing = TRUE), ]

# Prepare heatmap data
tissue_color <- read_excel('~/samplesheets/2021/20210222_Color_Scheme.xlsx') %>% 
  with(setNames(HEX, TISSUE))
meta <- read_excel("~/samplesheets/2021/20210226_MouseArray_SampleTableV4.xlsx")

# Remove extra columns and set row names
rownames(s) <- s$Sample_ID
s <- s[s$TISSUE != 'NA', ]
fs <- s[, 2:(ncol(s) - 3)]

# Rearrange so rows are samples
fs <- data.frame(t(fs), check.names = FALSE)

# Create heatmap
png("~/gallery/20210421_Igf2_ManualSorting_H19.png", height = 20, width = 15, units = 'in', res = 300)
visualizeGene('Igf2', fs, refversion = "mm10", platform = "MM285", 
              show.sampleNames = FALSE, dwstream = 80000, cluster.samples = FALSE) +
  WColorBarV(meta$Tissue_Corrected[match(colnames(fs), meta$IDAT)], 
             cmp = CMPar(label2color = tissue_color), LeftOf('betas')) +
  WLegendV(NULL, BottomRightOf('betas', h.pad = 0.1), height = .06)
dev.off()
