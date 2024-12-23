library(seqinr)
library(Rsamtools)

## Step 1: Reading the manifest probe CSV file
manifest_path <- '/mnt/isilon/zhoulab/labprojects/20200228_Mouse_Array_Project/20200721_manifestV2/mm10_LEGX_B3.manifest.sesame-base.cpg-sorted.csv.gz'
d <- read.csv(manifest_path)

## Step 2: Ensuring sequences are in the correct format
d$AlleleA_Probe_Sequence <- gsub('r', 'A', d$AlleleA_Probe_Sequence)
d$AlleleA_Probe_Sequence <- toupper(d$AlleleA_Probe_Sequence)

d$AlleleB_Probe_Sequence <- gsub('r', 'A', d$AlleleB_Probe_Sequence)
d$AlleleB_Probe_Sequence <- toupper(d$AlleleB_Probe_Sequence)

## Step 3: Writing FASTA files for sequences
seqA <- as.list(d$AlleleA_Probe_Sequence)
seqB <- as.list(d$AlleleB_Probe_Sequence)
id <- d$Probe_ID

fasta_output_path_A <- '/mnt/isilon/zhoulab/labprojects/20200228_Mouse_Array_Project/20200721_manifestV2/mm10_LEGX_B3.manifest.sesame-base.cpg-sorted_SequenceA.fasta'
fasta_output_path_B <- '/mnt/isilon/zhoulab/labprojects/20200228_Mouse_Array_Project/20200721_manifestV2/mm10_LEGX_B3.manifest.sesame-base.cpg-sorted_SequenceB.fasta'

write.fasta(seqA, id, fasta_output_path_A)
write.fasta(seqB, id, fasta_output_path_B)

## Linux commands for biscuit alignment (run outside R):
# cd /mnt/isilon/zhoulab/labprojects/20200228_Mouse_Array_Project/20200721_manifestV2/
# biscuit align -t 4 ~/references/mm10/biscuit/mm10.fa mm10_LEGX_B3.manifest.sesame-base.cpg-sorted_SequenceA | samtools view -O BAM -o mm10_LEGX_B3.manifest.sesame-base.cpg-sorted_MappedA.bam
# biscuit align -t 4 ~/references/mm10/biscuit/mm10.fa mm10_LEGX_B3.manifest.sesame-base.cpg-sorted_SequenceB | samtools view -O BAM -o mm10_LEGX_B3.manifest.sesame-base.cpg-sorted_MappedB.bam

## Step 4: Reading the mapped BAM files back into R
bamA <- data.frame(scanBam("/mnt/isilon/zhoulab/labprojects/20200228_Mouse_Array_Project/20200721_manifestV2/mm10_LEGX_B3.manifest.sesame-base.cpg-sorted_MappedA.bam"))
bamB <- data.frame(scanBam("/mnt/isilon/zhoulab/labprojects/20200228_Mouse_Array_Project/20200721_manifestV2/mm10_LEGX_B3.manifest.sesame-base.cpg-sorted_MappedB.bam"))

## Step 5: Renaming the columns
names(bamA) <- paste("AlleleA_MappedInfo", names(bamA), sep = "_")
names(bamB) <- paste("AlleleB_MappedInfo", names(bamB), sep = "_")

## Step 6: Combining the manifest data with mapped BAM data
df <- cbind(d, bamA, bamB)

## Step 7: Checking for unmapped probes based on design type
# Type I probes (flag 4 means unmapped)
t1 <- subset(df, DESIGN == "I" & AlleleA_MappedInfo_flag == 4)
cat("Unmapped Type I probes: ", nrow(t1), "\n")

# Type II probes can map to either Allele A or B, so check both for unmapped
t2 <- subset(df, DESIGN == "II" & AlleleA_MappedInfo_flag == 4 & AlleleB_MappedInfo_flag == 4)
cat("Unmapped Type II probes: ", nrow(t2), "\n")

## Step 8: Saving the results
saveRDS(df, file = "/mnt/isilon/zhoulab/labprojects/20200228_Mouse_Array_Project/20200721_manifestV2/mm10_LEGX_B3.manifest.sesame-base.cpg-sorted_Mapped.rds")

## Step 9: CpG Annotation Step
# Create a reference BED file extending TSS interval by 1500bp upstream
# Linux commands (run outside R):
# awk '$4=="+"{print $1,$2-1500,$3,$0;} $4=="-"{print $1,$2,$3+1500,$0}' /mnt/isilon/zhou_lab/projects/20191221_references/mm10/annotation/gtf/Mus_musculus.GRCm38.99.gtf.sorted.bed >> a.text
# awk '$2<0{$2=0}1' a.text >> CleanedReference.bed

## Step 10: Inside R: Process the cleaned reference file
x <- read.table("CleanedReference.bed")
x <- cbind(x[1:3], x[7:length(x)])  # Removing original columns

write.table(x, 'TSSeditedReference.bed', quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

## Step 11: Creating mapped bed file
mappedBed <- cbind(df[14], df[16], df[16] + 2, df$Probe_ID)
mappedBed <- na.omit(mappedBed)

write.table(mappedBed, file = "/mnt/isilon/zhoulab/labprojects/20200228_Mouse_Array_Project/20200721_manifestV2/mappedBed.bed", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

## Linux commands for bedtools intersection (run outside R):
# bedtools intersect -a TSSeditedReference.bed -b mappedBed.bed -wo -bed >> Reference_Matched_Probes.bed
# bedtools intersect -a Reference_Matched_Probes.bed -b ~/references/mm10/annotation/cpg_island/mm10_cpgIsland_FromUCSC.bed -wo -bed >> Reference+CPG_Matched_Probes.bed

## Step 12: Calculating distance to TSS (run outside R):
# awk '{if($4=="+"){print $2-$15;}if($4=="-"){print $16-$3;}}' Reference+CPG_Matched_Probes.bed >> TSS_distance.bed

## Step 13: Reading and combining annotation data
bed <- as.data.frame(read.table("/mnt/isilon/zhoulab/labprojects/20200228_Mouse_Array_Project/20200721_manifestV2/Reference+CPG_Matched_Probes.bed", header = FALSE, sep = "\t", stringsAsFactors = FALSE))
bed2 <- as.data.frame(read.table("/mnt/isilon/zhoulab/labprojects/20200228_Mouse_Array_Project/20200721_manifestV2/TSS_distance.bed", header = FALSE, sep = "\t", stringsAsFactors = FALSE))

df2 <- cbind(bed, bed2)

## Step 14: Renaming columns for clarity
names(df2) <- c('Transcript_Chr', 'Transcript_Start', 'Trnscript_End', 'transcript_chr',
                'transcriptoriginal_start', 'transcriptoriginal_end', 'Strand_Orientation', 
                'Annotation_Source', 'Annotation_Type', 'Gene_Name', 'Gene_Annotation', 
                'Annotation_Source', 'Ensemble_ID', 'Probe_Chr', 'Probe_Start', 'Probe_End', 
                'Probe_Name', 'Probe_Overlap', 'CPG_Chromosome', 'CPG_Start', 'CPG_End', 
                'CPG_Name', 'CPG_Overlap', 'TSS_Distance')

## Step 15: Saving the annotated data
saveRDS(df2, file = "/mnt/isilon/zhoulab/labprojects/20200228_Mouse_Array_Project/20200721_manifestV2/mm10_LEGX_B3.manifest.sesame-base.cpg-sorted_CPGandGenesAnnotated.rds")
