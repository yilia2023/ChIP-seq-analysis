# load necessary packages
library(tidyverse)
library(DESeq2)
library(magrittr)

# import data
con1 <- read.table("05_coverage/H3K27me3_Ctrl_rep1_count.tsv", sep = "\t",
                   col.names = c("chr", "start", "end", "con1", "a", "b", "c"))
con2 <- read.table("05_coverage/H3K27me3_Ctrl_rep2_count.tsv", sep = "\t",
                   col.names = c("chr", "start", "end", "con2", "a", "b", "c"))
LIF1 <- read.table("05_coverage/H3K27me3_noLIF_rep1_count.tsv", sep = "\t",
                   col.names = c("chr", "start", "end", "LIF1", "a", "b", "c"))
LIF2 <- read.table("05_coverage/H3K27me3_noLIF_rep2_count.tsv", sep = "\t",
                   col.names = c("chr", "start", "end", "LIF2", "a", "b", "c"))
con1 <- con1[, 1:4]
con2 <- con2[, 1:4]
LIF1 <- LIF1[, 1:4]
LIF2 <- LIF2[, 1:4]

# Prepare the data in the required format for DESeq2 analysis
con_peak <- merge(con1, con2, by=c("chr", "start", "end"), all = TRUE)
LIF_peak <- merge(LIF1, LIF2, by=c("chr", "start", "end"), all = TRUE)
count_peak <- merge(con_peak, LIF_peak, by=c("chr", "start", "end"), all = TRUE)
row.names(count_peak) <- paste(count_peak$chr, ":", count_peak$start, "-", count_peak$end)
countData = count_peak[,c("con1", "con2", "LIF1", "LIF2")]
row.names(countData) <- row.names(count_peak)

# creat colData
condition <- factor(c(rep("con", 2), rep("LIF", 2)), levels = c("con", "LIF"))
colData <- data.frame(row.names = colnames(countData), condition)

# creat DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ condition)

# differential peaks analysis by DEseq2
dds1 <- DESeq(dds)
res <- results(dds1)

# Sort by log2 Fold Change and convert to a data frame
resOrdered1 <- res[order(res$log2FoldChange),] %>% as.data.frame()

# identify differential peaks 
up <- resOrdered1[which(resOrdered1$log2FoldChange > 1 & resOrdered1$padj < 0.05),]    
down <- resOrdered1[which(resOrdered1$log2FoldChange < -1 & resOrdered1$padj < 0.05),]  

# Extract the chromatin position information
up_peak <- merge(count_peak[, 1:3], up, by = "row.names", all = FALSE)
down_peak <- merge(count_peak[, 1:3], down, by = "row.names", all = FALSE)

# export the results
write.table(up_peak[,2:4], file = "H3K27me3_peaks_down.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = F)
write.table(down_peak[,2:4], file = "H3K27me3_peaks_up.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = F)
