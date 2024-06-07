setwd("/home/vant/ATAC_seq/mouse_ATAC/peak_calling/htseq")

if (!requireNamespace("BiocManager", quietly = TRUE))
  
  install.packages("BiocManager")



BiocManager::install("DESeq2")

BiocManager::install("apeglm")

BiocManager::install("vsn",force = TRUE)

BiocManager::install("hexbin")

BiocManager::install("pheatmap")

library("apeglm")

library("DESeq2")

library("ggplot2")

library("vsn")


setwd("/home/vant/ATAC_seq/OSCC/peak_calling/htseq/")
directory = getwd()





##"greps" everything in current directory that has "htseq" in it for use as your files

#this methods assumes all htseq files are still separate and not together in a single file

sampleFiles <- data.frame(c("s8_O68_1_EKDL230002753-1A_sort_dedup_htseq_counts","S20_ATAC_O37_1_sort_dedup_htseq_counts","S26_ATAC_O43_NT_sort_dedup_htseq_counts","S27_ATAC_O43_NT_sort_dedup_htseq_counts","S28_ATAC_O43_T_sort_dedup_htseq_counts","S29_ATAC_O43_T_sort_dedup_htseq_counts","s3_O57_EKDL230002748-1A_HWGLWDSX5_sort_dedup_htseq_counts","s4_O64_EKDL230002749-1A_HWGLWDSX5_sort_dedup_htseq_counts","s7_O67_EKDL230002752-1A_HWGLWDSX5_sort_dedup_htseq_counts","s1_O55_1_EKDL230002746-1A_sort_dedup_htseq_counts","s2_O55_2_EKDL230002747-1A_sort_dedup_htseq_counts","S30_ATAC_O41_1_sort_dedup_htseq_counts","S31_ATAC_O41_2_sort_dedup_htseq_counts",
                            "s4_ATAC_O11_2_sort_dedup_htseq_counts","s5_O66_1_EKDL230002750-1A_sort_dedup_htseq_counts","s6_O66_2_EKDL230002751-1A_sort_dedup_htseq_counts","S1_ATAC_O18_1_sort_dedup_htseq_counts","S2_ATAC_O18_2_sort_dedup_htseq_counts","S5_ATAC_O22_1_sort_dedup_htseq_counts","S6_ATAC_O22_2_sort_dedup_htseq_counts","S11_ATAC_O27_1_sort_dedup_htseq_counts","S12_ATAC_O27_2_sort_dedup_htseq_counts","S7_ATAC_O24_1_sort_dedup_htseq_counts","S8_ATAC_O24_2_sort_dedup_htseq_counts","S19_ATAC_O35_1_sort_dedup_htseq_counts","s8_ATAC_O14_2_sort_dedup_htseq_counts","s9_ATAC_O14_3_sort_dedup_htseq_counts"))

colData <- data.frame(condition=factor(c("H","H","H","H","H","H","H","H","H","H","H","H","H","L","L","L","L","L","L","L","L","L","L","L","L","L","L")))


##specifies what condition the files are



sampleCondition <- colData

sampleTable <- data.frame(sampleName = sampleFiles, fileName = sampleFiles, condition = sampleCondition)



ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, design = ~ condition)

keep <- rowSums(counts(ddsHTSeq)) >= 1

ddsHTSeq <- ddsHTSeq[keep,]

ddsHTSeq$condition



dds <- DESeq(ddsHTSeq)

#actual DE analysis



#plotCounts(dds, intgroup = "condition", normalized = TRUE, transform = TRUE, xlab = "group")



dds <- estimateSizeFactors(dds)

normalizedcounts <- counts(dds, normalized=TRUE)

write.csv(as.data.frame(normalizedcounts), file = "Resultados/new_high_and_low/Normalized_DESeq_high_vs_low.csv")



res <- results(dds, contrast=c("condition","H","L"))

res

#tells DESeq2 what parameter to compare on and specifically which things in that parameter (in this case condition: A vs D)



write.csv(as.data.frame(res), file = "Resultados/new_high_and_low/OSCC_high_vs_low.csv")




resultsNames(dds)

resLFC <- lfcShrink(dds, coef="condition_L_vs_H", type="apeglm")

resLFC

resLFCOrdered <- resLFC[order(resLFC$pvalue),]

resLFCOrdered

#shrinking data (log2 transformation?)



resOrdered <- res[order(res$pvalue),]

resOrdered

res05 <- results(dds, alpha=0.05)


#reorder results by pvalue



summary(res)

sum(res$padj < 0.1, na.rm=TRUE) #na.rm means should NaN's be removed?

#see how many adjusted pvalues were less than 0.1; by default FDR is also 0.1 (make sure to match pvalue/FDR if changing)



res05 <- results(dds, alpha = 0.05)

summary(res05)

sum(res05$padj < 0.05, na.rm=TRUE)



plotMA(res, ylim=c(-2,2))

plotMA(resLFC, ylim=c(-2,2))

#plotting MA plots for either results or shrunken results (shrunken better for MA since it gets rid of 0 count noise on the left side)



idx <- identify(res$baseMean, res$log2FoldChange)

rownames(res)[idx]

#identify points on MA plot by clicking - start script, click, then finish/quit(esc) to see



plotCounts(dds,which.min(res$padj), intgroup = "condition")

#plot counts of gene with lowest adjusted pvalue in the "interesting group" condition



d <- plotCounts(dds, gene = which.min(res$padj), intgroup="condition",returnData=TRUE)

ggplot(d, aes(x=condition, y=count)) + geom_point(position=position_jitter(w=0.1, h=0)) + scale_y_log10(breaks=c(25))

#customize the plot, returnData returns only dataframe (no plot) for use in plotting with ggplot



mcols(res)$description

resSig_pval005 <- subset(resOrdered, pvalue < 0.05)

resSig_padj005 <- subset(resOrdered, padj < 0.05)

resSig


write.csv(as.data.frame(resOrdered), file = "Resultados/new_high_and_low/Deseq_OSCC_high_low_Ordered.csv")
write.csv(as.data.frame(resSig_pval005), file = "Resultados/new_high_and_low/Deseq_OSCC_high_low_pval005.csv")
write.csv(as.data.frame(resSig_padj005), file = "Resultados/new_high_and_low/Deseq_OSCC_high_low_padj005.csv")

#write results

#makes a subset table with just the things with padj < 0.1



vsd <- vst(dds, blind=FALSE)

rld <- rlog(dds, blind=FALSE)

#either variance stabilizing (parametric fit) or regularized log transformation of your data



ntd <- normTransform(dds)



meanSdPlot(assay(ntd))

meanSdPlot(assay(vsd))

meanSdPlot(assay(rld))

#look at how transformations affect mean SD's





library("pheatmap")


select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]

df <- as.data.frame(colData(dds)["condition"]) #can also include [,c("condition","type)] for multi parameter comparison

pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,cluster_cols=FALSE, annotation_col=df)

test<-assay(ntd)[select,]

pheatmap(test, cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE,, annotation_col=df)

pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         
         cluster_cols=FALSE, annotation_col=df)



plotPCA(vsd, intgroup=c("condition"))

pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  
  geom_point(size=3)


my_results_threshold <- results(object = dds,
                                contrast = c("condition", "H", "L"),
                          
                                alpha = 0.05,
                                pAdjustMethod = "BH",
                                tidy = TRUE
)
mat <- assay(vsd)[head(order(my_results_threshold$padj), 10), ] 
pheatmap(mat,cluster_rows=FALSE,cluster_cols=FALSE, annotation_col=df)

