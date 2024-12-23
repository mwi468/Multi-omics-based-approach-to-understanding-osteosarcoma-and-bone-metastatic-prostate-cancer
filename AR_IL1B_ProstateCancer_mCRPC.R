rm(list=ls())

library(ggplot2)
library(matrixStats)
library(survival)
library(data.table)
library(stringr)
library(pheatmap)
library(RColorBrewer)
library(tidyverse)

# location of gene TPM file 

x <- read.table("C:/Users/Waleed/Dropbox/SharedDesktopFiles/Lundberg_at_al/counts_and_TPMs/Counts_TPM_Genes.txt", header =TRUE, check.names = FALSE)

x_s <- x[2:(ncol(x)-1)]
x_s <- x_s[,names(x_s) %in% s]
x_s$gene <- x$gene_name

## gene names are in last column 

#NEPC
#genes <- c('SYP', 'CHGA', 'ENO2', 'RTN4RL2')

#original AR activity
genes <- c('KLK3','KLK2', 'FKBP5', 'STEAP1', 'STEAP2', 'PLPP1', 'RAB3B', 'NKX3-1', 'ACSL3','IL1B')

# New AR
#genes <- c('KLK3','NKX3-1', 'CHRNA2','SLC45A3','TRGC1','TRGC1','NAP1L2','IL1B')

#combined AR
#genes <- c('KLK3','NKX3-1', 'CHRNA2','SLC45A3','TRGC1','TRGC1','NAP1L2', 'KLK2', 'FKBP5', 'STEAP1', 'STEAP2', 'PLPP1', 'RAB3B','ACSL3','IL1B')

#genes2 <- c('KLK3','AR','IL1B')

#genes <- read.table("C:/Users/Waleed/Documents/ARv7_upregulated.txt")

x2 <- x_s[x_s$gene %in% genes,]
#exp_log <- log(x2[,1:(ncol(x2)-1)]+1)

exp2 <- data.frame(t(x2[1:(ncol(x2)-1)]))
#names(exp2) <- x2$gene_name

names(exp2) <- x2$gene

#exp2$mean <- rowMeans(as.matrix(exp2[,1:ncol(exp2)]))

#cor2 <- cor.test(exp2$RTN4RL2,exp2$mean, method= "spearman" )
#cor2

#ggplot(exp2,aes(mean,RTN4RL2)) + geom_point() + geom_smooth(method="lm") + ylab("NGR2 Log(TPM + 1)") + xlab("Mean NEPC_genes(SYP,CHGA,ENO2) Log(TPM + 1)") + ggtitle("Correlation of NEPC Activity vs NGR2 (n = 147)") + annotate("text", label = paste("Correlation:",round(cor2$estimate, digits=1),"Pvalue:",round(cor2$p.value, digits=5)), x = 3, y = .5) # + facet_zoom(ylim = c(0, 8))

#select AR expresssion

##
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

##

zexp <- data.frame(apply(exp2,2,Zscore));
rownames(zexp) <- rownames(exp2)

zexp$samples <- rownames(zexp)
zexp <- zexp %>% relocate(IL1B, .after = last_col())
zexp$mean <- rowMedians(as.matrix(zexp[,1:(ncol(zexp)-2)]))

cor2 <- cor.test(zexp$IL1B,zexp$mean, method= "spearman" )
cor2

ggplot(zexp,aes(IL1B,mean)) + geom_point() + geom_smooth(method="lm") + ylab("New Gene AR activity (z score)") + xlab("IL1B (zscore)") + ggtitle("AR Activity of combined gene set vs IL1B") + annotate("text", label = paste("Correlation:",round(cor2$estimate, digits=4),"Pvalue:",round(cor2$p.value, digits=9)), x = 2, y = 2) # + facet_zoom(ylim = c(0, 8))


IL-1β methylation analysis in mCRPC 
rm(list=ls())

library(data.table)
library(stringr)
library(pheatmap)
library(RColorBrewer)
library(tidyverse)
library(gridExtra)

## load data

# USCS gene expression data 
exp <- read.table("C:/Users/Waleed/Dropbox/SharedDesktopFiles/IL1B_collaboration_data/featurecounts_tpm_gene_gencodev28_n100.tsv", header = TRUE)

# USCS methylation data 
meth <- read.table("C:/Users/Waleed/Dropbox/SharedDesktopFiles/IL1B_collaboration_data/methylation_il1b_region.txt", header = TRUE)

# keep same samples 

exp.data <- exp[,2:ncol(exp)]
exp.names <- exp[,1]

meth.data <- meth[,3:ncol(meth)]
meth.position <- meth[,2]

exp <- exp.data[names(exp.data) %in% names(meth.data)]

meth <- meth.data[names(meth.data) %in% names(exp.data)]

# ordering meth tp be same as exp sample order (cols)

reorder_idx <- match(names(exp),names(meth))  
meth <- meth[,reorder_idx]

# add back feature id and cpg Position
exp <- cbind(exp.names,exp)
names(exp)[1] <- 'feature_id'
meth <- cbind(meth.position,meth)

## Subset genes 

# NEPC gene list from PMID: 30655535
#list <- read.table("C:/Users/u/Documents/NEPC genes.txt", sep = " ")

#AR activity genes 
argenes <- c('IL1B','AR','KLK3','NKX3-1', 'CHRNA2','SLC45A3','TRGC1','TRGC1','NAP1L2', 'KLK2', 'FKBP5', 'STEAP1', 'STEAP2', 'PLPP1', 'RAB3B','ACSL3')


# Getting Ensmble IDs from gencode V28
# Gencode version file was downloaded from https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/gencode.v28.annotation.gtf.gz
# this file was simplified using awk in shell  
# zcat gencode.v28.annotation.gtf.gz | awk 'BEGIN { OFS=FS="\t" } $3=="gene"{if($1=="MT")$1="M"; match($9, /gene_id "([^;]*) /, gene_id); match($9, /gene_name "([^;]*) /, gene_name); {print "AllGenes",gene_id[1], gene_name[1]}}'

names <- fread("C:/Users/Waleed/Dropbox/SharedDesktopFiles/gencode.v28.simplified.txt", header = FALSE)

## subset expression for relevant genes 
#subset names 
names2keep2 <- names[names$V3 %in% argenes,]
names(names2keep2)[2] <- 'feature_id'

## subset exp based on ensembl id in names2keep

exp2 <- exp[exp$feature_id %in% names2keep2$feature_id,]
exp2 <-  merge(names2keep2, exp, by = 'feature_id' )

exp <- exp2

# split column separating gene id and names 

n.exp <- exp[,4:ncol(exp)]
rownames(n.exp) <- exp$V3

## z score normalization


# for each Gene Z score normalize expression across samples 
# rows are genes, columns are samples

### convert 0 to really small number 
n.exp[n.exp == 0] = .0001

x <- as.vector(unlist(n.exp[9,]))
ZscoreRow = function(x,output){
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

zexp <- apply(n.exp,1,ZscoreRow);
zexp <- data.frame(zexp)

rownames(zexp) <- colnames(n.exp);
colnames(zexp) <- rownames(n.exp);

n.exp2 <- data.frame(t(zexp))
colnames(n.exp2) <- row.names(n.exp)

zexp <- zexp %>% relocate(IL1B, .after = last_col())
zexp <- zexp %>% relocate(AR, .after = last_col())

zexp$median <- rowMeans(as.matrix(zexp[,1:(ncol(zexp)-2)]))

n.exp3 <- data.frame(zexp[order(zexp$IL1B),])
rownames(n.exp3) <- rownames(zexp)
 
n.exp4 <- n.exp3[n.exp3$AR < 0,]
n.exp4 <- n.exp4[n.exp4$median < 0,]
#n.exp4 <- n.exp4[n.exp4$IL1B < 0,]

n.exp5 <- n.exp3[n.exp3$AR >= 0,]
n.exp5 <- n.exp5[n.exp5$median >= 0,]
#n.exp5 <- n.exp5[n.exp5$IL1B >= 0,]

# order meth data in same way 
options(scipen=999)

#cg20157753 (chr2:112,836,766-112,836,767), cg07935264(chr2:112,836,798-112,836,799), cg18773937 (chr2:112,837,034-112,837,035), cg23149881 (chr2:112,837,077-112,837,078), 

start1 <- 112836766 - 500
end1   <- 112836767 + 500

start2 <- 112836798 - 500
end2   <- 112836799 + 500

start3 <- 112837034 - 500
end3   <- 112837035 + 500

start4 <- 112837077 - 500
end4 <- 112837078 + 500

meth1 <- meth[meth$meth.position >= start1,]
meth1 <- meth1[meth1$meth.position <= end1,]

meth2 <- meth[meth$meth.position >= start2,]
meth2 <- meth2[meth2$meth.position <= end2,]

meth3 <- meth[meth$meth.position >= start3,]
meth3 <- meth3[meth3$meth.position <= end3,]

meth4 <- meth[meth$meth.position >= start4,]
meth4 <- meth4[meth4$meth.position <= end4,]

methc <- rbind(meth3,meth4)
methc <- distinct(methc)
#meth <- meth[meth$meth.position <= end,]

methf <- data.frame(t(methc))
n.e <- n.exp5
reorder_idx <- match(rownames(n.e),rownames(methf))
methf <- methf[reorder_idx,]
#names(n.e) <- ''
P1 <- pheatmap(n.e,scale = "none", colorRampPalette(rev(brewer.pal(n = 6, name = "RdYlBu")))(5), breaks = c(-1,-0.5,0,.05,1), cluster_cols = FALSE, cluster_rows = FALSE, cellheight = 4, cellwidth = 8, cex = 1 ,show_colnames= TRUE, show_rownames= FALSE, treeheight_row = 0, treeheight_col = 0, fontsize = 7, legend = TRUE, height = 10, width = 15, border_color=NA)
P2 <- pheatmap(methf, scale = "none", cluster_rows = FALSE, cluster_cols = FALSE, cellheight = 4, cellwidth = 1, cex = 1, show_colnames= FALSE, show_rownames= FALSE, treeheight_row = 0, treeheight_col = 0, fontsize = 7, legend = TRUE, height = 10, width = 5)
plot_list=list()
plot_list[['GeneExp']]=P1[[4]]
plot_list[['Meth_hm']]=P2[[4]]
P <- grid.arrange(grobs=plot_list, ncol=2, nrow=1)




ggsave(paste("NEPCmethylation_NEW.pdf"),P, width = 30, height = 12, dpi=300, scale=1)

cor.test(select$NEPCActivity, select$IL1B, method= "spearman" )

ggplot(select,aes(IL1B,NEPCActivity)) + geom_point() + geom_smooth(method="lm") + ylab("NEPC activity ") + xlab("NRG normalized expression") + ggtitle("Correlation of NEPC activity and NRG2 normalized expressions in Low AR, NEPC cases") + annotate("text", label = paste("Correlation:",round(cor$estimate, digits=3),"Pvalue:",round(cor$p.value, digits = 3)), x = 10, y = 10)
