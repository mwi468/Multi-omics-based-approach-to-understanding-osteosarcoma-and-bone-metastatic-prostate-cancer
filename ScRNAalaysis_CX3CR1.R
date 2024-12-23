library(data.table)
library(dplyr)

t <- fread("D:/PMID_30566875/GSE181294_scRNAseq.ano.csv.gz")
t2 <- fread("D:/PMID_30566875/GSE181294_slide.seq.ano.csv.gz")

## combined tumor file 

m <- fread("D:/PMID_30566875/Tumor_CX3CR1/combined.tsv")
#m <- m[ , colSums(is.na(m))==0]
names <- colnames(m)
m <- data.frame(m)
#m <- m[ , colSums(is.na(m))==0]
#m <- m %>% select_if(~ !any(is.na(.))) # remove NA
m <- data.frame(t(m))
colnames(m)[1] <- 'Values'
m$Names <- paste(unlist(names))
m <- m[ rowSums(is.na(m))==0,]

t <- data.frame(t)
t <- t[2:nrow(t),]

colnames(t) <- c('Names', 'Cell', 'Sample')

c <- merge(m,t, by = "Names")

## calculate Percentage of cells that express CX3CR1
c2 <- c[c$Values > 0,]

#test <- c[c$Cell == "Tumor",]

## combined normal file 
t <- data.frame(t)
t <- t[2:nrow(t),]

colnames(t) <- c('Names', 'Cell', 'Sample')

c <- merge(m,t, by = "Names")
test <- c[c$Cell == "Tumor",]


for (i in unique(c$Sample)) {
  test  = c[c$Sample == "SCG-PCA5-T-LG",]
  print(paste(i,nrow(all[all$Cell == "Tumor",])))
  }

table <- data.table(table(c2$Cell))
table2 <- data.table(table(c$Cell))

PercentCRexp <- data.table(cbind(table$V1,table$N/table2$N*100))
PercentCRexp$V2 <- as.numeric(PercentCRexp$V2)
PercentCRexp <- PercentCRexp[order(PercentCRexp$V2,decreasing = TRUE),]

pdf(file = "D:/Figures/PC_SCe_PMID30566875_CX3CR1percentage",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 5) # The height of the plot in inches

barplot(PercentCRexp$V2, main="Percentage with Cx3CR1 expression",names.arg = PercentCRexp$V1, las=2, cex.names=.75)

dev.off()

