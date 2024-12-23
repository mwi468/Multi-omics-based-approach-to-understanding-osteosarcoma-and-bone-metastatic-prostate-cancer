#clear current workspace
rm(list = ls())

library(readxl) # for reading xlsx files
library(Rtsne) # for doing tSNE dimensionality reduxtion
library(gridExtra) # for plorring in combined figure
library(tidyverse) # for ggplot as well as otheruseful R tools

#load the sample list
# we will sublist samples cleared for public sharing

All_Samples <- read_excel("/Users/iqbalw/Dropbox/ZhouLab/Lab Data/SampleSheets/2021/20210713_MouseArray_SampleTableFinal.xlsx")

#only want validated samples that will be shared online
#white list
R_Samples <- All_Samples[All_Samples$Data_publicly_shared == 1,]

# load the mouse array betas file
# samples are columns, probes are rows

betass <- readRDS("/Users/iqbalw/Documents/20210913_1740_mouse_array_betas.rds" )

#load the color dataframe for coloring the tSNE
color_df <- readRDS( "/Users/iqbalw/Downloads/20210914_mouse_array_samples_GroupColorCodes.rds")

#Subset both files for the same samples as in white list

WL_samples <- as.list(R_Samples$IDAT)

betass <- betass[,colnames(betass) %in% WL_samples]

color_df <- color_df[color_df$Sample_ID %in% WL_samples, ]

#incase diff samples or repeats make sure colordf and betass has the same samples

mat_keep_rows <- as.list(color_df$Sample_ID)

betass <- betass[,colnames(betass) %in% mat_keep_rows]

mat_keep_rows2 <- as.list(colnames(betass))

color_df <- color_df[color_df$Sample_ID %in% mat_keep_rows2, ]

# Order the samples in both to be in same order

betass <- betass[, order(colnames(betass)) ]

color_df <- color_df[order(color_df$Sample_ID),]


## Print number of samples and Probes used

print(paste("number of samples is: ", ncol(betass)))
print(paste("number of probes is: ", nrow(betass)))

# printing samples by category 

print(table(color_df$TISSUE))

print(table(color_df$EXP_Group))

print(table(color_df$STRAIN))

print(table(color_df$SEX))

#For tSNE features need to be columns, in this case probes
# make samples as rows and probes as columns

betass <- t(betass)

#Before filtering 
dim(betass)
#1119 293199

## Remove rows with more than 50% NA
betass <- betass[which(rowMeans(!is.na(betass)) > 0.5), ]

betass <- betass[, which(colMeans(!is.na(betass)) > 0.5)]

#After filtering 
dim(betass)
#1091 287304

# For other cases replace NA with col Means (probe average beta value score)

k <- which(is.na(betass), arr.ind=TRUE)
betass[k] <- colMeans(betass, na.rm=TRUE)[k[,1]]

#Only Want Cg Probes

betass <- betass[,colnames(betass) %like% "cg"]

dim(betass)

#1091 280588

t <- Rtsne(betass, dims=2, perplexity=75)

df <- data.frame(x = t$Y[,1], #remain the same across diff. graphs below
                 y = t$Y[,2]) #remain the same across diff. graphs below

## using GridExtra for plotting
#coloring the samples based on annotations and colors in color_df file

#going to use i %% 2 because sample annotations are in even columns, only want to focus on even columns
# i + 1 (in odd columns) is the corresponding color to use in the plot

#storing graphs for all the diff categories then plotting them togeter

#subset the color_df because after filtering some of the samples were removed
# error: replacement has 1119 rows, data has 1091



color_df <- color_df[color_df$Sample_ID %in% rownames(betass), ]


for (i in 1:ncol(color_df)){
  
  #tissue plot
  
  if((i %% 2) == 0 && names(color_df)[i] == 'TISSUE'){
    print('OK')
    df$colour <- unlist(color_df[i + 1]) #change this to diff cololor scheme: example Exp_merged2$Hex
    df$label <- unlist(color_df[i]) # change this to the grouping you want to color, example Exp_merged2$TISSUE
    df$label[is.na(df$label)] <- 'NA'
    p1<-ggplot(df, aes(x, y, colour = label)) + geom_point(size = .1) + scale_colour_manual(values=setNames(df$colour,df$label)) + labs(title=paste("tSNE_",names(color_df)[i]), x ="tSNE1", y = "tSNE2")  + theme(legend.position = "none")
    #ggsave(paste("/Users/iqbalw/tSNE/updated/tSNE_5",names(color_df)[i],".pdf"), width = 7.5, height = 5)
    
  }
  
  #cell line plot 
  
  else if((i %% 2) == 0 && names(color_df)[i] == 'Cell_Line') {
    print('OK')
    df$colour <- unlist(color_df[i + 1]) #change this to diff cololor scheme: example Exp_merged2$Hex
    df$label <- unlist(color_df[i]) # change this to the grouping you want to color, example Exp_merged2$TISSUE
    df$label[is.na(df$label)] <- 'NA'
    p2 <- ggplot(df, aes(x, y, colour = label)) + geom_point(size = .1) + scale_colour_manual(values=setNames(df$colour,df$label)) + labs(title=paste("tSNE_",names(color_df)[i]), x ="tSNE1", y = "tSNE2")  + theme(legend.position = "none")
    #ggsave(paste("/Users/iqbalw/tSNE/updated/tSNE_5",names(color_df)[i],".pdf"), width = 7.5, height = 5)
    
    
  }
  
  # experiment group plot
  
  else if((i %% 2) == 0 && names(color_df)[i] == 'EXP_Group') {
    print('OK')
    df$colour <- unlist(color_df[i + 1]) #change this to diff cololor scheme: example Exp_merged2$Hex
    df$label <- unlist(color_df[i]) # change this to the grouping you want to color, example Exp_merged2$TISSUE
    df$label[is.na(df$label)] <- 'NA'
    p3 <- ggplot(df, aes(x, y, colour = label)) + geom_point(size = .1) + scale_colour_manual(values=setNames(df$colour,df$label)) + labs(title=paste("tSNE_",names(color_df)[i]), x ="tSNE1", y = "tSNE2") + theme(legend.position = "none")
    #ggsave(paste("/Users/iqbalw/tSNE/updated/tSNE_5",names(color_df)[i],".pdf"), width = 7.5, height = 5)
    
    
  }
  
  # STRAIN plot 
  
  else if((i %% 2) == 0 && names(color_df)[i] == 'STRAIN') {
    print('OK')
    df$colour <- unlist(color_df[i + 1]) #change this to diff cololor scheme: example Exp_merged2$Hex
    df$label <- unlist(color_df[i]) # change this to the grouping you want to color, example Exp_merged2$TISSUE
    df$label[is.na(df$label)] <- 'NA'
    p4 <- ggplot(df, aes(x, y, colour = label)) + geom_point(size = .1) + scale_colour_manual(values=setNames(df$colour,df$label)) + labs(title=paste("tSNE_",names(color_df)[i]), x ="tSNE1", y = "tSNE2") + theme(legend.position = "none")
    #ggsave(paste("/Users/iqbalw/tSNE/updated/tSNE_5",names(color_df)[i],".pdf"), width = 7.5, height = 5)
    
    
  }
  #sex plot 
  
  else if((i %% 2) == 0 && names(color_df)[i] == 'SEX') {
    print('OK')
    df$colour <- unlist(color_df[i + 1]) #change this to diff cololor scheme: example Exp_merged2$Hex
    df$label <- unlist(color_df[i]) # change this to the grouping you want to color, example Exp_merged2$TISSUE
    df$label[is.na(df$label)] <- 'NA'
    p5 <- ggplot(df, aes(x, y, colour = label)) + geom_point(size = .1) + scale_colour_manual(values=setNames(df$colour,df$label)) + labs(title=paste("tSNE_",names(color_df)[i]), x ="tSNE1", y = "tSNE2") + theme(legend.position = "none")
    #ggsave(paste("/Users/iqbalw/tSNE/updated/tSNE_5",names(color_df)[i],".pdf"), width = 7.5, height = 5)
    
    
  }
  # Age plot
  
  else if((i %% 2) == 0 && names(color_df)[i] == 'Age') {
    print('OK')
    df$colour <- unlist(color_df[i + 1]) #change this to diff cololor scheme: example Exp_merged2$Hex
    df$label <- unlist(color_df[i]) # change this to the grouping you want to color, example Exp_merged2$TISSUE
    df$label[is.na(df$label)] <- 'NA'
    p6 <- ggplot(df, aes(x, y, colour = label)) + geom_point(size = .1) + scale_colour_manual(values=setNames(df$colour,df$label)) + labs(title=paste("tSNE_",names(color_df)[i]), x ="tSNE1", y = "tSNE2") + theme(legend.position = "none")
    #ggsave(paste("/Users/iqbalw/tSNE/updated/tSNE_5",names(color_df)[i],".pdf"), width = 7.5, height = 5)
    
  }
  #Mean Meth Plot
  
  else if((i %% 2) == 0 && names(color_df)[i] == 'meanMeth') {
    print(i)
    df$colour <- unlist(color_df[i + 1]) #change this to diff cololor scheme: example Exp_merged2$Hex
    df$label <- unlist(color_df[i]) # change this to the grouping you want to color, example Exp_merged2$TISSUE
    df$label[is.na(df$label)] <- 'NA'
    p7 <- ggplot(df, aes(x, y, colour = label)) + geom_point(size = .1) + scale_colour_manual(values=setNames(df$colour,df$label)) + labs(title=paste("tSNE_",names(color_df)[i]), x ="tSNE1", y = "tSNE2") + theme(legend.position = "none")
    #ggsave(paste("/Users/iqbalw/tSNE/updated/tSNE_5",names(color_df)[i],".pdf"), width = 7.5, height = 5)
    
  }
  #PCG meth plot
  
  else if((i %% 2) == 0 && names(color_df)[i] == 'pcgMeth') {
    print(i)
    df$colour <- unlist(color_df[i + 1]) #change this to diff cololor scheme: example Exp_merged2$Hex
    df$label <- unlist(color_df[i]) # change this to the grouping you want to color, example Exp_merged2$TISSUE
    df$label[is.na(df$label)] <- 'NA'
    p8 <- ggplot(df, aes(x, y, colour = label)) + geom_point(size = .1) + scale_colour_manual(values=setNames(df$colour,df$label)) + labs(title=paste("tSNE_",names(color_df)[i]), x ="tSNE1", y = "tSNE2") + theme(legend.position = "none")
    #ggsave(paste("/Users/iqbalw/tSNE/updated/tSNE_5",names(color_df)[i],".pdf"), width = 7.5, height = 5)
    
  }
  # Tumor Vs Normal Plot
  
  else if((i %% 2) == 0 && names(color_df)[i] == 'TumorVsNormal') {
    print(i)
    df$colour <- unlist(color_df[i + 1]) #change this to diff cololor scheme: example Exp_merged2$Hex
    df$label <- unlist(color_df[i]) # change this to the grouping you want to color, example Exp_merged2$TISSUE
    df$label[is.na(df$label)] <- 'NA'
    p9 <- ggplot(df, aes(x, y, colour = label)) + geom_point(size = .1) + scale_colour_manual(values=setNames(df$colour,df$label)) + labs(title=paste("tSNE_",names(color_df)[i]), x ="tSNE1", y = "tSNE2") + theme(legend.position = "none")
    #ggsave(paste("/Users/iqbalw/tSNE/updated/tSNE_5",names(color_df)[i],".pdf"), width = 7.5, height = 5)
    
  }
  
  else {
  }
}

## Print the graphs into a combined file at desired location

#tissue (Figure 4A left), followed by sex (Figure 4A right), strain, and age

#P1 tissue 
#P2 Cell line
#p3 Exp Group
#p4 strain
#p5 sex
#p6 age
#p7 mean Meth
#p8 PCG meth
#P9 TvN

P <- grid.arrange(p1,p5,p4,p6,p2,p9,p3,p7,p8,nrow = 2, ncol= 5, widths = c(2, 2, 2, 2, 2))

ggsave(paste("~/gallery/20210811_WhiteList_tSNE_All.pdf"), P, width = 7, height = 3, dpi=300, scale=1.5)

