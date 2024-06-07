data(sdc,package='BulkSignalR')
normal <- grep("^N", names(sdc))
bsrdm <- prepareDataset(sdc[,-normal])

# define the comparison

bsrdm.comp <- as.BSRDataModelComp(bsrdm)
colA <- as.integer(1:5)
colB <- as.integer(8:15)

# As an example here, we generate random values 
# but user should provide his own logFC and
# associated pvalues from DGE ouputs. 

n <- nrow(ncounts(bsrdm.comp))
stats <- data.frame(pval=runif(n), logFC=rnorm(n, 0, 2),
                    expr=runif(n, 0, 10))
rownames(stats) <- rownames(ncounts(bsrdm.comp))


# we first define the cluster comparison and add it
# to the BSRDataModelComp object.

bsrcc <- defineClusterComp(bsrdm.comp, colA, colB, stats)
bsrdm.comp <- addClusterComp(bsrdm.comp, bsrcc, "random.example")

# finally we infer ligand-receptor interactions from the comparison

bsrinf.comp <- initialInference(bsrdm.comp, max.pval=1,"random.example")

bsrinf.redP <- reduceToPathway(bsrinf.comp)
bsrinf.redBP    <- reduceToBestPathway(bsrinf.comp)
bsrinf.L    <- reduceToLigand(bsrinf.comp)
bsrinf.R    <- reduceToReceptor(bsrinf.comp)
bsrinf.redP  <- reduceToPathway(bsrinf.comp)  
bsrinf.redPBP <- reduceToBestPathway(bsrinf.redP)



test<-caf30 %>% filter(row.names(caf30) %in% names) 
