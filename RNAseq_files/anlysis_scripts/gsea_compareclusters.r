#############################################
# GSEA Analysis with compareClusters        #
#                                           #
# - Author: Álvaro Martínez Martínez        #
#############################################

# First, load all packages and objects from DE Analysis
library(edgeR)
library(Biobase)
library(PCAtools)
library(NOISeq)
library(knitr)
library(dplyr)
# These packages include E. coli db
library(clusterProfiler)
library(enrichplot)
library(GOSemSim)
library(org.EcK12.eg.db)

load(file="C13_noiseq.rda")
load(file="TT_alive.rda")
load(file="TT_dead.rda")


##########################
# BUILD GO SEMSIM MATRIX #
##########################

# Executing this function takes some time, so we can directly load
# the semantic similarity matrix once computed
#go.sem <- godata("org.EcK12.eg.db", ont = "BP")
#save(go.sem, file="GOSem.rda")
load("GOSem.rda")

####################
# GSEA ALIVE CAS13 #
####################
# Extract topTable object from the list
# Codes: 1 = Lwa, 2 = Rfx, 3 = Lbu, 4 = eLbu
tt1.new <- tt.list[[1]]
tt2.new <- tt.list[[2]]
tt3.new <- tt.list[[3]]
tt4.new <- tt.list[[4]]
# Arrange a sorted vector based on -log10(raw.p) signed by the logFC.
# Repeat this for each Cas13 condition
stat1 <- ifelse(tt1.new$table$logFC < 0, log10(tt1.new$table$PValue), -log10(tt1.new$table$PValue))
names(stat1) <- rownames(tt1.new)
stat1 <- sort(stat1, decreasing = TRUE)

stat2 <- ifelse(tt2.new$table$logFC < 0, log10(tt2.new$table$PValue), -log10(tt2.new$table$PValue))
names(stat2) <- rownames(tt2.new)
stat2 <- sort(stat2, decreasing = TRUE)

stat3 <- ifelse(tt3.new$table$logFC < 0, log10(tt3.new$table$PValue), -log10(tt3.new$table$PValue))
names(stat3) <- rownames(tt3.new)
stat3 <- sort(stat3, decreasing = TRUE)

stat4 <- ifelse(tt4.new$table$logFC < 0, log10(tt4.new$table$PValue), -log10(tt4.new$table$PValue))
names(stat4) <- rownames(tt4.new)
stat4 <- sort(stat4, decreasing = TRUE)

# Create a list with the sorted vectors and 
# execute compareCluster with GSEA (gseGO)
alive <- list(Lwa = stat1, Rfx = stat2, Lbu = stat3, eLbu = stat4)
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# THIS FUNCTION TAKES A WHILE TO EXECUTE
# IT'S BEST TO SAVE THE RESULT AFTERWARDS
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#alive.res <- compareCluster(alive, fun = "gseGO",
#                            ont = "BP",
#                            OrgDb = org.EcK12.eg.db,
#                            keyType = 'SYMBOL',
#                            verbose = T,
#                            minGSSize = 10,
#                            maxGSSize = 500,
#                            pvalueCutoff = .05,
#                            pAdjustMethod = "BH",
#                            seed = T)
#save(alive.res, file = "CC_alive.rda")
load("CC_alive.rda")

# With this you can reduce redundancy based on semsim data:
alive.res <- simplify(alive.res, cutoff = 0.5, semData = go.sem)

# Some basic plots of the GO terms enriched
dotplot(alive.res, showCategory = 50)
treeplot(pairwise_termsim(alive.res, semData = go.sem), showCategory = 100, offset = rel(2),
         offset_tiplab = rel(2), nWords = 3, hclust_method = "average")


###################
# GSEA DEAD CAS13 #
###################
# Extract topTable object from the list
# Codes: 1 = Lwa, 2 = Rfx, 3 = Lbu, 4 = eLbu
tt1.dead <- dead.tt.list[[1]]
tt2.dead <- dead.tt.list[[2]]
tt3.dead <- dead.tt.list[[3]]
tt4.dead <- dead.tt.list[[4]]
# Arrange a sorted vector based on -log10(raw.p) signed by the logFC.
# Repeat this for each Cas13 condition
stat1.d <- ifelse(tt1.dead$table$logFC < 0, log10(tt1.dead$table$PValue), -log10(tt1.dead$table$PValue))
names(stat1.d) <- rownames(tt1.dead)
stat1.d <- sort(stat1.d, decreasing = TRUE)

stat2.d <- ifelse(tt2.dead$table$logFC < 0, log10(tt2.dead$table$PValue), -log10(tt2.dead$table$PValue))
names(stat2.d) <- rownames(tt2.dead)
stat2.d <- sort(stat2.d, decreasing = TRUE)

stat3.d <- ifelse(tt3.dead$table$logFC < 0, log10(tt3.dead$table$PValue), -log10(tt3.dead$table$PValue))
names(stat3.d) <- rownames(tt3.dead)
stat3.d <- sort(stat3.d, decreasing = TRUE)

stat4.d <- ifelse(tt4.dead$table$logFC < 0, log10(tt4.dead$table$PValue), -log10(tt4.dead$table$PValue))
names(stat4.d) <- rownames(tt4.dead)
stat4.d <- sort(stat4.d, decreasing = TRUE)

# Create a list with the sorted vectors and 
# execute compareCluster with GSEA (gseGO)
dead <- list(dLwa = stat1.d, dRfx = stat2.d, dLbu = stat3.d, deLbu = stat4.d)
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# THIS FUNCTION TAKES A WHILE TO EXECUTE
# IT'S BEST TO SAVE THE RESULT AFTERWARDS
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#dead.res <- compareCluster(dead, fun = "gseGO",
#                           ont = "BP",
#                           OrgDb = org.EcK12.eg.db,
#                           keyType = 'SYMBOL',
#                           verbose = T,
#                           minGSSize = 10,
#                           maxGSSize = 500,
#                           pvalueCutoff = .05,
#                           pAdjustMethod = "BH",
#                           seed = T)
#save(dead.res, file = "CC_dead.rda")
load("CC_dead.rda")

# With this you can reduce redundancy based on semsim data:
dead.res <- simplify(dead.res, cutoff = 0.5, semData = go.sem)

# Some basic plots of the GO terms enriched
dotplot(dead.res, showCategory = 50)
treeplot(pairwise_termsim(dead.res, semData = go.sem), showCategory = 100, offset = rel(2),
         offset_tiplab = rel(4), nWords = 3, hclust_method = "average")
