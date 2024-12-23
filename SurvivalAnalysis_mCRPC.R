rm(list=ls())

library(ggplot2)
library(matrixStats)
library(survival)
library(survminer)
library(data.table)
library(dplyr)
library(reshape)

######
survival <- read.table("C:/Users/Waleed/Dropbox/SharedDesktopFiles/Lundberg_at_al//wcdt_fu_data_231128.txt", sep = " ", header = TRUE)
#samples to keep 

#########################################
#### load needed data 
x <- read.table("C:/Users/Waleed/Dropbox/SharedDesktopFiles/Lundberg_at_al//site_wcdt.csv", sep =",", header = TRUE)

xa <- read.table("C:/Users/Waleed/Dropbox/SharedDesktopFiles/Lundberg_at_al/counts_and_TPMs/Counts_TPM_Genes.txt", header =TRUE, check.names = FALSE)
xb <- read.table("C:/Users/Waleed/Dropbox/SharedDesktopFiles/Lundberg_at_al/counts_and_TPMs/Counts_TPM_Genes_PROMOTE.txt", header =TRUE, check.names = FALSE)

x1 <- merge(xa[1:(ncol(xa)-1)],xb[1:(ncol(xb))],by="FEATURE_ID")
rm(xa,xb)

#stemness we chose to look at median z score of SALL4, NANOG, OCT4a (POU5F1) and ALDH1A1 (ALDH1)
### subset data for needed genes
genes <- c('KLK3','KLK2', 'FKBP5', 'STEAP1', 'STEAP2', 'PLPP1', 'RAB3B', 'NKX3-1', 'ACSL3','IL1B')
#genes <- c('KLK3','NKX3-1', 'CHRNA2','SLC45A3','TRGC1','TRGC1','NAP1L2', 'KLK2', 'FKBP5', 'STEAP1', 'STEAP2', 'PLPP1', 'RAB3B','ACSL3', 'IL1B')
#genes <- c('SYP', 'CHGA', 'ENO2','IL1B')
#genes <- read.table("C:/Users/Waleed/Dropbox/SharedDesktopFiles/NEPCGenesUP_PMID29132337.csv", sep = ",")
#genes <- read.table("C:/Users/Waleed/Dropbox/SharedDesktopFiles/NEPC genes.txt", sep = " ")
#genes <- genes$V1
#genes[89] <- 'IL1B'

x2 <- x1[x1$gene_name %in% unlist(genes),]


##

##########################################
############## Z score



n.exp <- x2[,2:(ncol(x2)-1)]
rownames(n.exp) <- x2$gene_name
n.exp <- data.frame(t(n.exp))

zeros_in_columns <- sapply(n.exp, function(x) sum(x == 0))
print(zeros_in_columns)

n.exp <- log2(n.exp+1)
#zexp <- data.frame(scale(n.exp, center = FALSE, scale = TRUE))
zexp <- data.frame(apply(n.exp, 2, scale))
rownames(zexp) <- rownames(n.exp)

## remove NEPC
#zexp$sampleID <- paste(rownames(zexp))
#zexp <- zexp %>%
#  mutate(PatientID = sub("^(\\w+-\\d+)-.*", "\\1", sampleID))

## 
#filtered_df <- survival %>% filter(Subtypes == "AR-/NE+" | Subtypes == "AR+/NE+")

## subset 
#x_NEPC <- zexp[zexp$PatientID %in% filtered_df$PATIENT_ID,]
#x_NEPCremoved <- zexp[!zexp$PatientID %in% filtered_df$PATIENT_ID,]

## bone
#x <- read.table("C:/Users/Waleed/Downloads/site_wcdt.csv", sep =",", header = TRUE)
#subset_x <- x[x$Site_detailed == 'Bone',]

#x_NEPCremoved <- x_NEPCremoved[x_NEPCremoved$PatientID %in% subset_x$PATIENT_ID,]

## AR activity

zexp <- zexp %>% relocate(IL1B, .after = last_col())


zexp$mean <- rowMeans(zexp[,1:(ncol(zexp)-1)])
zexp$median <- rowMedians(as.matrix(zexp[,1:(ncol(zexp)-2)]))


# Load necessary library
library(ggplot2)

# Define thresholds
lower_threshold <-  median(combined$IL1B)
upper_threshold <- median(combined$IL1B +.0001)

# Create z_data from zexp$mean if not already done
z_data <- data.frame(z_scores = combined$IL1B)

# Categorize data based on thresholds
z_data$category <- cut(
  z_data$z_scores,
  breaks = c(-Inf, lower_threshold, upper_threshold, Inf),
  labels = c("Low", "Normal", "High"),
  right = FALSE
)

# Calculate counts in each category
category_counts <- table(z_data$category)

# Create labels with counts for the legend
category_labels <- paste(
  names(category_counts),
  " (n=",
  as.numeric(category_counts),
  ")",
  sep = ""
)

# Create a named vector for colors
category_colors <- c("Low" = "red", "Normal" = "gray", "High" = "blue")

# Create density data frame
density_data <- density(z_data$z_scores)
density_df <- data.frame(x = density_data$x, y = density_data$y)

# Categorize the x-values in the density data
density_df$category <- cut(
  density_df$x,
  breaks = c(-Inf, lower_threshold, upper_threshold, Inf),
  labels = c("Low", "Normal", "High"),
  right = FALSE
)

# Calculate positions for annotations
annotation_positions <- data.frame(
  category = names(category_counts),
  x = c(
    mean(density_df$x[density_df$category == "Low"]),
    mean(density_df$x[density_df$category == "Normal"]),
    mean(density_df$x[density_df$category == "High"])
  ),
  y = max(density_df$y) * 0.8,
  label = paste(names(category_counts), "\n(n=", as.numeric(category_counts), ")", sep = "")
)

# Plot the density curve with counts in legend and annotations
ggplot(density_df, aes(x = x, y = y, fill = category)) +
  geom_area(alpha = 0.5) +
  geom_text(
    data = annotation_positions,
    aes(x = x, y = y, label = label),
    color = "black",
    size = 5,
    hjust = 0.5
  ) +
  scale_fill_manual(
    values = category_colors,
    name = "Category",
    labels = category_labels
  ) +
  labs(
    title = "Density Plot of Mean Z-Scores with Counts",
    x = "Z-Score",
    y = "Density"
  ) +
  theme_minimal()

# Low AR activity

hi <- zexp[zexp$mean > 1,]
#low <- zexp[zexp$IL1B > 0,]
#hi <- hi[hi$IL1B < median(combined$IL1B),] # low IL1B
hi <- hi[hi$IL1B < median(combined$IL1B),] # low IL1B
hi$ARactivity <- 'High'

low <- zexp[zexp$mean < median(combined$mean),]
#low <- zexp[zexp$IL1B > 0,]
#low <- low[low$IL1B > median(combined$IL1B),] # high IL1B
low <- low[low$IL1B > median(combined$IL1B),] # high IL1B
low$ARactivity <- 'Low'

combined <-rbind(hi,low)

ggplot(combined, aes(x=ARactivity, y=IL1B)) + 
  geom_violin(trim=TRUE)

## survival

# Load necessary libraries
library(survival)
library(survminer)
library(ggplot2)

# Subset the survival data based on 'hi' and 'low' groups
survivalh <- survival[survival$SAMPLE_ID %in% rownames(hi), ]
survivalh$ARactivity <- 'High'

survivall <- survival[survival$SAMPLE_ID %in% rownames(low), ]
survivall$ARactivity <- 'Low'

# Combine the two groups
survival2 <- rbind(survivalh, survivall)

# Fit the survival curves
fit <- survfit(Surv(time_from_biopsy_months, event) ~ ARactivity, data = survival2)

# Fit the Cox proportional hazards model
cox_model <- coxph(Surv(time_from_biopsy_months, event) ~ ARactivity, data = survival2)

# Summarize the Cox model
s <- summary(cox_model)
print(s)

# Customize the survival plot
surv_plot <- ggsurvplot(
  fit,
  data = survival2,
  pval = TRUE,                   # Add p-value from log-rank test
  conf.int = TRUE,               # Add confidence intervals
  risk.table = TRUE,             # Include risk table
  risk.table.col = "strata",     # Color risk table by strata
  palette = c("blue", "red"),    # Set custom colors
  xlab = "Time (Months)",        # X-axis label
  ylab = "Survival Probability", # Y-axis label
  legend.title = "AR Activity",  # Legend title
  legend.labs = c("High", "Low"),# Legend labels
  surv.median.line = "hv",       # Add median survival lines
  ggtheme = theme_minimal()      # Use a minimal theme
)

# Further customize the plot using ggplot2 functions
surv_plot$plot <- surv_plot$plot +
  theme(
    text = element_text(family = "Arial"),   # Change font
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.position = "bottom"               # Position legend
  ) +
  scale_color_manual(
    values = c("High" = "blue", "Low" = "red"),
    labels = c("High AR Activity", "Low AR Activity")
  )

# Customize the risk table
surv_plot$table <- surv_plot$table +
  theme(
    text = element_text(family = "Arial"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

# Print the customized plot
print(surv_plot)

