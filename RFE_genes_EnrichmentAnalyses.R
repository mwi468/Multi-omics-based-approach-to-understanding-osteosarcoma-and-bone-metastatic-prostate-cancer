rm(list=ls())

library(msigdbr)
library(fgsea)
library(clusterProfiler)
library(org.Hs.eg.db)

## load merged dataset and datasets

merged_df <- readRDS("C:/Users/u/Dropbox/SharedDesktopFiles/RFE_geneComparisonTable.rds")
hallmark_genesets <- msigdbr(species = "Homo sapiens", category = "H")
geneSets <- split(hallmark_genesets$entrez_gene, hallmark_genesets$gs_name)

############# Upregulated IDS

sub <- merged_df[merged_df$MestastaticBone == "Up",]
sub <- sub[order(sub$PrimaryBone, decreasing=TRUE), ]
sub$rank <- as.numeric(sub$rank)
sub_ranked <- sub[order(sub$rank), ]

## load 
## convert IDS
#### gene table with gene names and RFE rank for GSEA analysis

subset_table <- sub_ranked %>%
  select(Gene, rank)

geneList <- setNames(subset_table$rank, subset_table$Gene)

# Sort the geneList in decreasing order (required for GSEA)
geneList <- sort(geneList, decreasing = TRUE)

gene_entrez <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = names(geneList),    # Your gene symbols (names of geneList)
  column = "ENTREZID",       # Map to Entrez IDs
  keytype = "SYMBOL",        # Input type: SYMBOL
  multiVals = "first"        # In case of duplicates, take the first match
)

# Remove NAs and prepare the gene list
geneList2 <- geneList[!is.na(gene_entrez)]
names(geneList2) <- gene_entrez[!is.na(gene_entrez)]

####### Pathway enrichment analyses 
## Gene Ontology 
## takes in gene Symbols (use subset_table$Gene)

go_results <- enrichGO(
  gene = subset_table$Gene,         # Your gene list
  OrgDb = org.Hs.eg.db,     # Organism database
  keyType = "SYMBOL",       # Input gene ID type
  ont = "BP",               # Ontology: BP (Biological Process), MF, or CC
  pAdjustMethod = "BH",     # Adjust p-values for multiple testing
  pvalueCutoff = 0.05,      # Significance threshold
  qvalueCutoff = 0.2        # Adjusted p-value threshold
)

head(go_results)

## KEGG 

kegg_results <- enrichKEGG(
  gene = names(geneList2),         # Your gene list
  organism = "hsa",         # Organism code for humans
  keyType = "kegg",         # Input gene ID type (SYMBOL/ENTREZID)
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)

# View results
head(kegg_results)

## MSigDB gene sets 

hallmark_results <- enricher(
  gene = names(geneList2),
  TERM2GENE = hallmark_genesets[, c("gs_name", "entrez_gene")],
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)

head(hallmark_results)

# make the dotplots 

dotplot_plot1 <- dotplot(go_results, showCategory = 10) +
  ggtitle("Gene Ontology Analysis in Upregulated Genes") +
  theme_minimal()

dotplot_plot2 <- dotplot(kegg_results, showCategory = 10) +
  ggtitle("KEGG Analysis in Upregulated Genes") +
  theme_minimal()

dotplot_plot3 <- dotplot(hallmark_results, showCategory = 10) +
  ggtitle("MSigDB Hallmark Analysis in Upregulated Genes") +
  theme_minimal()

# Save the plot as a PDF with specified dimensions
ggsave("C:/Users/u/Dropbox/SharedDesktopFiles/GO_Up_dotplot.pdf", plot = dotplot_plot1, width = 10, height = 5, device = "pdf")
ggsave("C:/Users/u/Dropbox/SharedDesktopFiles/KEGG_Up_dotplot.pdf", plot = dotplot_plot2, width = 10, height = 5, device = "pdf")
ggsave("C:/Users/u/Dropbox/SharedDesktopFiles/MSigDB_Up_dotplot.pdf", plot = dotplot_plot3, width = 10, height = 5, device = "pdf")

##### downregulated genes 

sub <- merged_df[merged_df$MestastaticBone == "Down",]
sub <- sub[order(sub$PrimaryBone, decreasing=TRUE), ]
sub$rank <- as.numeric(sub$rank)
sub_ranked <- sub[order(sub$rank), ]

subset_table <- sub_ranked %>%
  select(Gene, rank)

geneList <- setNames(subset_table$rank, subset_table$Gene)

# Sort the geneList in decreasing order (required for GSEA)
geneList <- sort(geneList, decreasing = TRUE)

gene_entrez <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = names(geneList),    # Your gene symbols (names of geneList)
  column = "ENTREZID",       # Map to Entrez IDs
  keytype = "SYMBOL",        # Input type: SYMBOL
  multiVals = "first"        # In case of duplicates, take the first match
)

# Remove NAs and prepare the gene list
geneList2 <- geneList[!is.na(gene_entrez)]
names(geneList2) <- gene_entrez[!is.na(gene_entrez)]

####### Pathway enrichment analyses 

## only GO Molecular Function (MF) had results
## Gene Ontology 
## takes in gene Symbols (use subset_table$Gene)

go_results <- enrichGO(
  gene = subset_table$Gene,         # Your gene list
  OrgDb = org.Hs.eg.db,     # Organism database
  keyType = "SYMBOL",       # Input gene ID type
  ont = "MF",               # Ontology: BP (Biological Process), MF, or CC
  pAdjustMethod = "BH",     # Adjust p-values for multiple testing
  pvalueCutoff = 0.05,      # Significance threshold
  qvalueCutoff = 0.2        # Adjusted p-value threshold
)
head(go_results)

## KEGG 
kegg_results <- enrichKEGG(
  gene = names(geneList2),         # Your gene list
  organism = "hsa",         # Organism code for humans
  keyType = "kegg",         # Input gene ID type (SYMBOL/ENTREZID)
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)
# View results
head(kegg_results)

## MSigDB gene sets 
hallmark_results <- enricher(
  gene = names(geneList2),
  TERM2GENE = hallmark_genesets[, c("gs_name", "entrez_gene")],
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)

head(hallmark_results)

# make the dotplots 
dotplot_plot4 <- dotplot(go_results, showCategory = 10) +
  ggtitle("Gene Ontology Analysis in Downregulated Genes") +
  theme_minimal()

## nothing found for the KEGG
dotplot_plot5 <- dotplot(kegg_results, showCategory = 10) +
  ggtitle("KEGG Analysis in Upregulated Genes") +
  theme_minimal()

## nothing foungdin MSIgDB
dotplot_plot6 <- dotplot(hallmark_results, showCategory = 10) +
  ggtitle("MSigDB Hallmark Analysis in Upregulated Genes") +
  theme_minimal()

ggsave("C:/Users/u/Dropbox/SharedDesktopFiles/GO_Down_dotplot.pdf", plot = dotplot_plot4, width = 10, height = 5, device = "pdf")

