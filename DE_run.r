############################################################
# Differential expression analysis
############################################################


library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(ggrepel)
library(reshape2)
library(lmerTest)
library(MAST)
library(data.table)

setwd("/mnt/e/Shared/lk/rerun/SS")
source("mast.de.functions.r")

########################################################

test.formula <- formula(~n_genes + tissue + procedure + Layer + Type + (1 | PubID) + (1 | PubIDSample))
test.formula.prefilter <- formula(~n_genes + tissue + procedure + Layer + Type)
test.contrasts <- c("TypeN", "TypeH")
subsample.pattern <- list(~Type, ~PubID, ~PubIDSample)
max.cells <- 10000
min_expression_frac <- 0.1

de.dir <- "de-results"
prefilter.P.thresh <- 0.05

#################################################################################
## Stromal ~5h

seur <- readRDS(file="scRNA_str.rds")

# Exclude inflamed cells
Idents(seur) <- seur@meta.data$status
seur <- subset(seur, idents=c("F","H","N"))

# modify tissue
seur@meta.data$tissue[seur@meta.data$tissue=="small intestine"] <- "small_intestine"
seur@meta.data$tissue[seur@meta.data$tissue=="ascending colon"] <- "ascending_colon"
table(seur@meta.data$tissue)

unique(seur@meta.data$annotation2)

out_label <- "Stromal"
for (celltype in unique(seur@meta.data$annotation2)) {
	run.celltype(seur, out_label, celltype)
}

#################################################################################
## Epithelial ~12h

seur <- readRDS(file="scRNA_epi.rds")

# Exclude inflamed cells
Idents(seur) <- seur@meta.data$status
seur <- subset(seur, idents=c("F","H","N"))

# modify tissue
seur@meta.data$tissue[seur@meta.data$tissue=="small intestine"] <- "small_intestine"
seur@meta.data$tissue[seur@meta.data$tissue=="ascending colon"] <- "ascending_colon"
table(seur@meta.data$tissue)

unique(seur@meta.data$annotation2)

out_label <- "Epithelial"
for (celltype in unique(seur@meta.data$annotation2)) {
	run.celltype(seur, out_label, celltype)
}


#################################################################################
## Immune ~12h

seur <- readRDS(file="scRNA_imm.rds")

# Exclude inflamed cells
Idents(seur) <- seur@meta.data$status
seur <- subset(seur, idents=c("F","H","N"))

# modify tissue
seur@meta.data$tissue[seur@meta.data$tissue=="small intestine"] <- "small_intestine"
seur@meta.data$tissue[seur@meta.data$tissue=="ascending colon"] <- "ascending_colon"
table(seur@meta.data$tissue)

unique(seur@meta.data$annotation2)

out_label <- "Immune"
for (celltype in unique(seur@meta.data$annotation2)) {
	run.celltype(seur, out_label, celltype)
}

########################################################
# after running all DE results
# put all DE results together

setwd(de.dir)
de.files <- list.files(pattern="*.csv")
length(de.files) # 65

de.cmb <- list()
for (fname in de.files) {
	fn <- file.path(de.dir, fname)
	if (file.exists(fn)) {
		de.res <- read.csv(fn)
		de.res$X <- substr(fname, 1, nchar(fname)-4)
		de.cmb <- c(de.cmb, list(de.res))
	}
}

de.cmb <- do.call(rbind, de.cmb)
de.cmb$loc <- gsub("^(\\w+)\\..*$", "\\1", de.cmb$X)
de.cmb$cy <- gsub("^\\w+\\.", "", de.cmb$X)
de.cmb$FDR.D <- p.adjust(ifelse(is.na(de.cmb$Pr..Chisq.D), de.cmb$Pr..Chisq.D.pre, de.cmb$Pr..Chisq.D), method='fdr')
de.cmb$FDR.C <- p.adjust(ifelse(is.na(de.cmb$Pr..Chisq.C), de.cmb$Pr..Chisq.C.pre, de.cmb$Pr..Chisq.C), method='fdr')
de.cmb$FDR <- p.adjust(ifelse(is.na(de.cmb$Pr..Chisq.), de.cmb$Pr..Chisq..pre, de.cmb$Pr..Chisq.), method='fdr')
de.cmb$coefD <- ifelse(is.na(de.cmb$coefD), de.cmb$coefD.pre, de.cmb$coefD)


de.cmb$coefD <- -1*de.cmb$coefD
de.cmb$direct.D <- de.cmb$coefD<0 # FALSE=UP, TRUE=DOWN 

# correct compartment
de.cmb$X[de.cmb$cy=="KCNN3−hi fibroblasts"] <- "Stromal.KCNN3−hi fibroblasts"
de.cmb$loc[de.cmb$cy=="KCNN3−hi fibroblasts"] <- "Stromal"
de.cmb$X[de.cmb$cy=="KRT-hi cells"] <- "Immune.KRT-hi cells"
de.cmb$loc[de.cmb$cy=="KRT-hi cells"] <- "Immune"

saveRDS(de.cmb, "de.cmb.rds")