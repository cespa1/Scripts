BiocManager::install("ggtree")
BiocManager::install("orthogene")
devtools::install_github("jcolinge/BulkSignalR")
install.packages("igraph")
install.packages("ragg")

library(BulkSignalR)
library(igraph)
library(dplyr)

##Cargar los pvalores y el logFC de las muestras
setwd("R-L/")
caf30_stats<-read.csv("Caf30_DEG_padj005.csv", header = TRUE, sep = ",")
names<-caf30_stats$X
caf30_stats<-read.csv("Caf30_DEG_padj005.csv", header = TRUE, row.names = names, sep = ",")
caf30_stats<-select(caf30_stats,-X)
caf30_stats<-data.frame(pval=caf30_stats$pvalue,logFC=caf30_stats$log2FoldChange)
row.names(caf30_stats)<-names

##Cargar los counts de los genes del RNAseq
caf30<-read.csv("Caf30.csv", header = TRUE, sep = ",")
names1<-caf30$X
caf30<-read.csv("Caf30.csv", header = TRUE, row.names = names1, sep = ",")
caf30<-select(caf30,-X)

bsrdm <- prepareDataset(caf30)
bsrdm.comp <- as.BSRDataModelComp(bsrdm)
colA<-c(1:3)
ColB<-c(4:6)
bsrcc <- defineClusterComp(bsrdm.comp, colA, colB, caf30_stats)


bsrdm <- learnParameters(bsrdm, 
                         plot.folder = "./", filename = "test")
bsrinf <- initialInference(bsrdm)

bsrinf.redP <- reduceToPathway(bsrinf)
bsrinf.redBP    <- reduceToBestPathway(bsrinf)
bsrinf.L    <- reduceToLigand(bsrinf)
bsrinf.R    <- reduceToReceptor(bsrinf)
bsrinf.redP  <- reduceToPathway(bsrinf)  
bsrinf.redPBP <- reduceToBestPathway(bsrinf.redP)

bsrsig.redBP <- getLRGeneSignatures(bsrinf.redBP, qval.thres=0.001)

scoresLR <- scoreLRGeneSignatures(bsrdm, bsrsig.redBP,
                                  name.by.pathway=FALSE)

simpleHeatmap(scoresLR[1:20,],
              path = "./",
              filename = "preliminar_data_scoresLR",
              column.names = TRUE,
              height = 20, width = 10,
              pointsize = 8,
              hcl.palette = "Cividis"                  
)

bsrsig.redPBP <- getLRGeneSignatures(bsrinf.redPBP, qval.thres=0.001)

scoresPathway <- scoreLRGeneSignatures(bsrdm, bsrsig.redPBP,
                                       name.by.pathway=TRUE)

simpleHeatmap(scoresPathway[1:50,],
              path = "./figures/",
              filename = "preliminar_data_scoresPathway",
              column.names = TRUE,
              width = 15,
              height = 10,
              pointsize = 12,
              hcl.palette = "Blue-Red 2"
)

pathway1 <- "Signaling by NOTCH4"
signatureHeatmaps(
  pathway = pathway1,
  bsrdm = bsrdm,
  bsrsig = bsrsig.redPBP,
  path = "./figures/",
  filename = "test_signatureheatmap",
  h.width  = 15,
  h.height = 10 ,
  show_column_names = TRUE)

alluvialPlot(bsrinf,
             keywords = c("COL4A1"),
             type = "L",
             qval.thres = 0.001,
             path = "./figures/",
             filename = "preliminar_data_alluvial",
             width  = 16,
             height = 12
)

#pathway1 <- "PD-1 signaling"
pathway2 <- "Interferon gamma signaling"
pathways <- c(pathway1,pathway2)

bubblePlotPathwaysLR(bsrinf,
                     pathways = pathways,
                     qval.thres  = 0.00001,
                     path = "./figures/",
                     color = "red",
                     filename  = "preliminar_data_bubble",
                     width  = 20,
                     height = 10,
                     pointsize = 8
                     #filter.L = c("ADAM12")
                     #filter.R = c("ITGA3")
)



chordDiagramLR (bsrinf,
                path = "./figures/",
                filename = "preliminar_data_chord",
                pw.id.filter = "R-HSA-202733",
                limit = 20,
                width = 5,
                height = 4.5
)