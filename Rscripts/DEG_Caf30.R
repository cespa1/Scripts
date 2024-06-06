install.packages("BiocManager")

BiocManager::install("DESeq2")

BiocManager::install("apeglm")

BiocManager::install("vsn",force = TRUE)

aBiocManager::install("hexbin")

BiocManager::install("pheatmap")

library("apeglm")
library("DESeq2")
library("ggplot2")
library("vsn")
library(dplyr)


setwd("R-L/")
caf30<-read.csv("Caf30.csv")
names<-caf30$X
caf30<-read.csv("Caf30.csv", row.names = names)
caf30<-select(caf30,-X)

nombre<-c("Cc","Cc","Cc","NoCc","NoCc","NoCc")
a<-colnames(caf30)
coldata<-data.frame(
  "tratamiento"=c("Cc","Cc","Cc","NoCc","NoCc","NoCc"),row.names = a
)

dds <- DESeqDataSetFromMatrix(countData = caf30,
                              colData = coldata,
                              design = ~ tratamiento)

dds <- DESeq(dds)
res <- results(dds)
res
write.csv(as.data.frame(res),file="Caf30_DEG.csv")
res1 <- results(dds)
res1

plotCounts(dds, gene=which.min(res$padj), intgroup="tratamiento")

ntd <- normTransform(dds)
##FALTA ORDENAR "RES" POR PVALOR Y PONER UN THRESHOLD EN 0.05


