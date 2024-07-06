#############################################
# DE Analysis with NOISeq Correction Script #
#                                           #
# - Author: Álvaro Martínez Martínez        #
#############################################

# Load libraries and ESet
library(edgeR)
library(Biobase)
library(PCAtools)
library(NOISeq)
library(knitr)
library(dplyr)

load(file = "../matrix_and_eset/C13_exprs.rda")

# Convert to dge
dge = DGEList(counts = exprs(C13.eSet), group = pData(C13.eSet)$condicion)

# Discard low count rows from it
thold <- 10/(min(dge$samples$lib.size)/1000000)
tokeep <- rowSums(cpm(exprs(C13.eSet)) > thold) >= 4
dge <- dge[tokeep,,keep.lib.size = FALSE]
# Discard rRNA genes (some got a little bit of counts, 
# prob due to bad clean-up during library prep)
rRNA <- grep("rr(s|l|f)", rownames(dge$counts))
dge = dge[-rRNA,,keep.lib.size = FALSE]
# Apply discard
C13.eSet <- C13.eSet[tokeep,]
C13.eSet <- C13.eSet[-rRNA,]

# ARSyNSeq Correction from NOISeq
C13.eSet.noise = ARSyNseq(C13.eSet, factor = "condicion", batch = FALSE, norm = "tmm", logtransf = FALSE)
# Save an ESet with NOISeq Correction applied
save(C13.eSet.noise, file="C13_noiseq.rda")

############
# OPTIONAL #
############
##################################################
# Correction can be compared in a PCA plot

# This is with data pre-corrected
prenorm.vals <- log(exprs(C13.eSet) + 1)
pca.res = PCAtools::pca(mat = prenorm.vals, metadata = dge$samples, scale = TRUE)
screeplot(pca.res)
biplot(pca.res, colby = "group", lab = dge$samples$group)
biplot(pca.res, showLoadings = TRUE, colby="group", lab = NULL, labSize = 5, pointSize = 5, sizeLoadingsNames = 5)
# This is with corrected data
norm.vals <- log(exprs(C13.eSet.noise) + 1)
pca.res = PCAtools::pca(mat = norm.vals, metadata = dge$samples, scale = TRUE)
screeplot(pca.res)
biplot(pca.res, colby = "group", lab = dge$samples$group)
biplot(pca.res, showLoadings = TRUE, colby="group", lab = NULL, labSize = 5, pointSize = 5, sizeLoadingsNames = 5)

# As it can be seen, noise is eliminated and samples of a same biological condition
# are more concentrated, but the altogether distribution does not substantially change
###################################################

##########################
# DE Analysis with edgeR #
##########################
dge = DGEList(counts = exprs(C13.eSet.noise), group = pData(C13.eSet.noise)$condicion)

###############
# ALIVE CAS13 #
###############
# First, run the analysis with samples concerning alive Cas13
new.dge <- dge[,grep("A(1|3|4|5|6|7|8|9|10)_",colnames(dge$counts))]
# Model matrix with the interesting samples
new.mm <- model.matrix(~ 0 + new.dge$samples$group)
colnames(new.mm) <- levels(new.dge$samples$group)
# Estimate commondisp and do the fit
new.dge <- estimateCommonDisp(new.dge)
fit.new <- glmFit(new.dge, design = new.mm)
# Define contrasts and make LRT
new.con = makeContrasts(conLWa = Lwa_gRNA - Lwa_empty, 
                        conRX = Rfx_gRNA - Rfx_empty,
                        conLbu = Lbu_gRNA - Lbu_empty,
                        coneLbu = eLbu_gRNA - eLbu_empty,
                        levels = new.mm)
# 1 = Lwa, 2 = Rfx, 3 = Lbu, 4 = eLbu
lrt1.new <- glmLRT(fit.new, contrast = new.con[,"conLWa"])
tt1.new <- topTags(lrt1.new, n=nrow(new.dge$counts))
sum(tt1.new$table$FDR <= .05)
lrt2.new <- glmLRT(fit.new, contrast = new.con[,"conRX"])
tt2.new <- topTags(lrt2.new, n=nrow(new.dge$counts))
sum(tt2.new$table$FDR <= .05)
lrt3.new <- glmLRT(fit.new, contrast = new.con[,"conLbu"])
tt3.new <- topTags(lrt3.new, n=nrow(new.dge$counts))
sum(tt3.new$table$FDR <= .05)
lrt4.new <- glmLRT(fit.new, contrast = new.con[,"coneLbu"])
tt4.new <- topTags(lrt4.new, n=nrow(new.dge$counts))
sum(tt4.new$table$FDR <= .05)

# Lots of DE genes. We can distinguish up and down regulated:
LWa.up <- which(tt1.new$table$FDR <= .05 & tt1.new$table$logFC > 0)
LWa.down <- which(tt1.new$table$FDR <= .05 & tt1.new$table$logFC < 0)
RX.up <- which(tt2.new$table$FDR <= .05 & tt2.new$table$logFC > 0)
RX.down <- which(tt2.new$table$FDR <= .05 & tt2.new$table$logFC < 0)
Lbu.up <- which(tt3.new$table$FDR <= .05 & tt3.new$table$logFC > 0)
Lbu.down <- which(tt3.new$table$FDR <= .05 & tt3.new$table$logFC < 0)
eLbu.up <- which(tt4.new$table$FDR <= .05 & tt4.new$table$logFC > 0)
eLbu.down <- which(tt4.new$table$FDR <= .05 & tt4.new$table$logFC < 0)

sumup <- data.frame(LWa=c(length(LWa.up), length(LWa.down)),
                    RX=c(length(RX.up), length(RX.down)),
                    Lbu=c(length(Lbu.up), length(Lbu.down)),
                    eLbu=c(length(eLbu.up), length(eLbu.down)),
                    row.names = c("Up-regulated", "Down-regulated")
)
kable(sumup)

# Finally, saving TopTables objects in a list
tt.list <- list(Lwa = tt1.new, 
                Rfx = tt2.new, 
                Lbu = tt3.new, 
                eLbu = tt4.new)
save(tt.list, file = "TT_alive.rda")

##############
# DEAD CAS13 #
##############
# First, run the analysis with samples concerning dead Cas13
dead.dge <- dge[,grep("A(1|11|12|13|14|15|16|17|18)_",colnames(dge$counts))]
# Model matrix with the interesting samples
dead.mm <- model.matrix(~ 0 + dead.dge$samples$group)
colnames(dead.mm) <- levels(dead.dge$samples$group)
# Estimate commondisp and do the fit
dead.dge <- estimateCommonDisp(dead.dge)
fit.dead <- glmFit(dead.dge, design = dead.mm)
# Define contrasts and make LRT
dead.con = makeContrasts(conLWa = dLwa_gRNA - dLwa_empty, 
                         conRX = dRfx_gRNA - dRfx_empty,
                         conLbu = dLbu_gRNA - dLbu_empty,
                         coneLbu = deLbu_gRNA - deLbu_empty,
                         levels = dead.mm)
# 1 = dLwa, 2 = dRfx, 3 = dLbu, 4 = deLbu
lrt1.dead <- glmLRT(fit.dead, contrast = dead.con[,"conLWa"])
tt1.dead <- topTags(lrt1.dead, n=nrow(dead.dge$counts))
sum(tt1.dead$table$FDR <= .05)
lrt2.dead <- glmLRT(fit.dead, contrast = dead.con[,"conRX"])
tt2.dead <- topTags(lrt2.dead, n=nrow(dead.dge$counts))
sum(tt2.dead$table$FDR <= .05)
lrt3.dead <- glmLRT(fit.dead, contrast = dead.con[,"conLbu"])
tt3.dead <- topTags(lrt3.dead, n=nrow(dead.dge$counts))
sum(tt3.dead$table$FDR <= .05)
lrt4.dead <- glmLRT(fit.dead, contrast = dead.con[,"coneLbu"])
tt4.dead <- topTags(lrt4.dead, n=nrow(dead.dge$counts))
sum(tt4.dead$table$FDR <= .05)

# Up/down genes
dLWa.up <- which(tt1.dead$table$FDR <= .05 & tt1.dead$table$logFC > 0)
dLWa.down <- which(tt1.dead$table$FDR <= .05 & tt1.dead$table$logFC < 0)
dRX.up <- which(tt2.dead$table$FDR <= .05 & tt2.dead$table$logFC > 0)
dRX.down <- which(tt2.dead$table$FDR <= .05 & tt2.dead$table$logFC < 0)
dLbu.up <- which(tt3.dead$table$FDR <= .05 & tt3.dead$table$logFC > 0)
dLbu.down <- which(tt3.dead$table$FDR <= .05 & tt3.dead$table$logFC < 0)
deLbu.up <- which(tt4.dead$table$FDR <= .05 & tt4.dead$table$logFC > 0)
deLbu.down <- which(tt4.dead$table$FDR <= .05 & tt4.dead$table$logFC < 0)

sumup.d <- data.frame(dLWa=c(length(dLWa.up), length(dLWa.down)),
                    dRX=c(length(dRX.up), length(dRX.down)),
                    dLbu=c(length(dLbu.up), length(dLbu.down)),
                    deLbu=c(length(deLbu.up), length(deLbu.down)),
                    row.names = c("Up-regulated", "Down-regulated")
)
kable(sumup.d)

# Finally, saving TopTables objects in a list
dead.tt.list <- list(dLwa = tt1.dead, 
                     dRfx = tt2.dead, 
                     dLbu = tt3.dead, 
                     deLbu = tt4.dead)
save(dead.tt.list, file = "TT_dead.rda")

