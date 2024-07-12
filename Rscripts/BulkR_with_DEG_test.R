caf30_stats.new <-caf30_stats %>% filter(row.names(caf30_stats) %in% test)
##Preparar el data set de counts 
bsrdm <- prepareDataset(caf30.new)

# define the comparison

bsrdm.comp <- as.BSRDataModelComp(bsrdm)
colA <- as.integer(1:3)
colB <- as.integer(4:6)

# we first define the cluster comparison and add it
# to the BSRDataModelComp object.

bsrcc <- defineClusterComp(bsrdm.comp, colA, colB, caf30_stats.new)
bsrdm.comp <- addClusterComp(bsrdm.comp, bsrcc, "random.example")

# finally we infer ligand-receptor interactions from the comparison

bsrinf <- initialInference(bsrdm.comp, max.pval=1,"random.example")

bsrinf.redP <- reduceToPathway(bsrinf.comp)
bsrinf.redBP    <- reduceToBestPathway(bsrinf.comp)
bsrinf.L    <- reduceToLigand(bsrinf.comp)
bsrinf.R    <- reduceToReceptor(bsrinf.comp)
bsrinf.redP  <- reduceToPathway(bsrinf.comp)  
bsrinf.redPBP <- reduceToBestPathway(bsrinf.redP)



test<-caf30 %>% filter(row.names(caf30) %in% names) 
