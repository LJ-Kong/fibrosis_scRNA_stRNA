###########################################################
# Raw data was preprocessed using CellRanger and Scanpy by 
# default settings. Fig. 1B, 1G and 1F were generated as 
# part of this
############################################################
# Figure 1C PCoAs
############################################################

library(ggplot2)
library(cowplot)
library(reshape2)
library(labdsv)

# Load metadata
df <- read.csv("./data/updated_SS_meta.csv")
dim(df) #352697     42
length(unique(df$annotation2)) # 77 

df <- df[grep("^Contaminant", df$annotation2, invert=T),] # exlude "Contaminant"
dim(df) #347017     42
length(unique(df$annotation2)) # 68 

propmat1 <- dcast(df, annotation2 ~ sampleid)

propmat <- matrix(as.numeric(unlist(propmat1[,-1])), nrow=nrow(propmat1)) # convert to numeric
rownames(propmat) <- propmat1$annotation2
colnames(propmat) <- colnames(propmat1)[-1]

propmat <- sweep(propmat, 2, colSums(propmat), FUN="/")
D <- dsvdis(t(propmat), index="bray/curtis")
pc <- pco(D)

dfc <- df[!duplicated(df$sampleid),]
dfc <- dfc[match(colnames(propmat),dfc$sampleid),]
dfc$pc1 <- pc$points[,1]
dfc$pc2 <- pc$points[,2]

varexp <- pc$eig / sum(pmax(0, pc$eig))

p1 <- ggplot(data=dfc) +
	geom_point(aes(x=pc1, y=pc2, fill=status), shape=21, size=3, color="black", stroke=0.2) + theme_bw() +
	scale_fill_manual(values=c(H="#95d962", N="#628cd9", F="#d9628c", I="#ffde05")) +
	xlab(sprintf("PCo 1 (%.1f%%)", 100*varexp[1])) +
	ylab(sprintf("PCo 2 (%.1f%%)", 100*varexp[2])) + ggtitle("Status") +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

p2 <- ggplot(data=dfc) +
	geom_point(aes(x=pc1, y=pc2, fill=fraction), shape=21, size=3, color="black", stroke=0.2) + theme_bw() +
	scale_fill_manual(values=c(epi="#42cbf5", imu="#f5428d")) +
	xlab(sprintf("PCo 1 (%.1f%%)", 100*varexp[1])) +
	ylab(sprintf("PCo 2 (%.1f%%)", 100*varexp[2])) + ggtitle("Fraction") +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

p3 <- ggplot(data=dfc) +
	geom_point(aes(x=pc1, y=pc2, fill=tissue), shape=21, size=3, color="black", stroke=0.2) + theme_bw() +
	scale_fill_manual(values=c("small intestine"="#ffde05", "colon"="#fc5d00", "ascending colon"="purple")) +
	xlab(sprintf("PCo 1 (%.1f%%)", 100*varexp[1])) +
	ylab(sprintf("PCo 2 (%.1f%%)", 100*varexp[2])) + ggtitle("Tissue location") +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

p4 <- ggplot(data=dfc) +
	geom_point(aes(x=pc1, y=pc2, fill=procedure), shape=21, size=3, color="black", stroke=0.2) + theme_bw() +
	scale_fill_manual(values=c("resection"="#5e9900", "biopsy"="#ffa1bb")) +
	xlab(sprintf("PCo 1 (%.1f%%)", 100*varexp[1])) +
	ylab(sprintf("PCo 2 (%.1f%%)", 100*varexp[2])) + ggtitle("Tissue type") +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

p <- plot_grid(p2, p3, p1, p4, nrow=1, rel_widths=c(3.9,4.5,3.8,4.1))
pdf("fig_Pco_20x4_horizon6.pdf", 20, 4)
print(p)
dev.off()

############################################################
# Figure 1D-E  Patient-wise PCoAs and
# Weighted Averages Scores for cell types
############################################################

library(ggplot2)
library(cowplot)
library(ggrepel)
library(reshape2)
library(labdsv)


# Load metadata
df <- read.csv("./data/updated_SS_meta.csv")
dim(df) #352697     42
length(unique(df$annotation2)) # 77

df <- df[grep("^Contaminant", df$annotation2, invert=T),] # exlude "Contaminant"
dim(df) #347017     42
length(unique(df$annotation2)) # 68 

##### keep patients that have all epi and imu
df2 <- read.csv("./data/combine_samples_by_fraction_manual.csv")
dim(df2) #102   6

# remove samples that are not paired
paired_samples <- df2$sampleid[df2$paired=="Y"]
df <- df[df$sampleid%in%paired_samples,]
df$patient_status <- paste(df$patient, df$status, sep="_")

propmat1 <- dcast(df, annotation2 ~ sampleid)
propmat <- matrix(as.numeric(unlist(propmat1[,-1])), nrow=nrow(propmat1)) # convert to numeric
rownames(propmat) <- propmat1$annotation2
colnames(propmat) <- colnames(propmat1)[-1]

propmat <- sweep(propmat, 2, colSums(propmat), FUN="/")

# add propotions together per patient
colnames(propmat) <- df2$patient_status[match(colnames(propmat), df2$sampleid)]
propmat <- sapply(split.default(as.data.frame(propmat), colnames(propmat)), rowSums, na.rm = TRUE)

D <- dsvdis(t(propmat), index="bray/curtis")
pc <- pco(D)

dfc <- df[!duplicated(df$patient_status),]
dfc <- dfc[match(colnames(propmat),dfc$patient_status),]
dfc$pc1 <- pc$points[,1]
dfc$pc2 <- pc$points[,2]

varexp <- pc$eig / sum(pmax(0, pc$eig))

p1 <- ggplot(data=dfc) +
	geom_point(aes(x=pc1, y=pc2, fill=status), shape=21, size=3, color="black", stroke=0.2) + theme_bw() +
	scale_fill_manual(values=c(H="#95d962", N="#628cd9", F="#d9628c", I="#ffde05")) +
	xlab(sprintf("PCo 1 (%.1f%%)", 100*varexp[1])) +
	ylab(sprintf("PCo 2 (%.1f%%)", 100*varexp[2])) + ggtitle("Status") +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

p3 <- ggplot(data=dfc) +
	geom_point(aes(x=pc1, y=pc2, fill=tissue), shape=21, size=3, color="black", stroke=0.2) + theme_bw() +
	scale_fill_manual(values=c("small intestine"="#ffde05", "colon"="#fc5d00", "ascending colon"="purple")) +
	xlab(sprintf("PCo 1 (%.1f%%)", 100*varexp[1])) +
	ylab(sprintf("PCo 2 (%.1f%%)", 100*varexp[2])) + ggtitle("Tissue location") +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

p4 <- ggplot(data=dfc) +
	geom_point(aes(x=pc1, y=pc2, fill=procedure), shape=21, size=3, color="black", stroke=0.2) + theme_bw() +
	scale_fill_manual(values=c("resection"="#5e9900", "biopsy"="#ffa1bb")) +
	xlab(sprintf("PCo 1 (%.1f%%)", 100*varexp[1])) +
	ylab(sprintf("PCo 2 (%.1f%%)", 100*varexp[2])) + ggtitle("Tissue type") +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

p <- plot_grid(p3, p1, p4, nrow=1, rel_widths=c(4.5,3.8,4.1))
pdf("fig_Pco_horizon_patientwise.pdf", 15, 4)
print(p)
dev.off()


#####################################################
# Weighted Averages Scores for cell types

library(vegan)
library(ggrepel)

wsc <- as.data.frame(wascores(pc$points, t(propmat)))
wsc$slope <- rowSums((wsc)^2)
wsc <- wsc[order(wsc$slope, decreasing=T),]
colnames(wsc) <- c("x", "y", "disxy")
wsc$annotation2 <- rownames(wsc)

p <- ggplot(data=dfc) + geom_point(aes(x=pc1, y=pc2, color=status), shape=16, size=4) + theme_bw() +
		scale_color_manual(values=c(H="#95d962", N="#628cd9", F="#d9628c", I="#ffde05")) +
		xlab(sprintf("PCo 1 (%.1f%%)", 100*varexp[1])) +
		ylab(sprintf("PCo 2 (%.1f%%)", 100*varexp[2])) + ggtitle("Status") +
		geom_text_repel(data=wsc[1:25,], aes(x=x, y=y, label=annotation2), size=4, box.padding = 0.5, color="black") + coord_equal() +
		geom_segment(data=wsc[1:25,], aes(x=0, y=0, xend=x, yend=y), arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="red")

pdf("fig_WeightedAveragesScoresTop25_patientwise.pdf", 12, 8)
print(p)
dev.off()



############################################################
# Figure 1F  Composion tests between resection vs. biopsy
############################################################

library(lmerTest)
library(MASS)
library(egg)
library(DirichletReg)
library(Matrix)
library(data.table)
library(tidyverse)
library(tidyr)
library("RColorBrewer")
library(cowplot)
library(scales)

source("dirichlet_functions.r")

# Load metadata
meta <- read.csv("./data/updated_SS_meta.csv")
dim(meta) #352697     42

# Exclude inflamed cells
meta <- meta[meta$status %in% c("F","H","N"),]
dim(meta) #345621     42

meta <- meta[grep("^Contaminant", meta$annotation2, invert=T),] # exlude "Contaminant"
dim(meta) #339990     42

meta$tissue[meta$tissue=="ascending colon"] <- "colon"
table(meta$tissue)

###########################
# create a new column change anno2 from factor to character
meta$anno22 <- as.character(meta$annotation2)

epi.freq = as.matrix(as.data.frame.matrix(table(meta$sampleid[meta$compart=="Epithelial cells"], meta$anno22[meta$compart=="Epithelial cells"])))
fib.freq = as.matrix(as.data.frame.matrix(table(meta$sampleid[meta$compart=="Stromal cells"], meta$anno22[meta$compart=="Stromal cells"])))
imm.freq = as.matrix(as.data.frame.matrix(table(meta$sampleid[meta$compart=="Immune cells"], meta$anno22[meta$compart=="Immune cells"])))

epi.tests <- cell_proportions_procedure(epi.freq, "imu")
fib.tests <- cell_proportions_procedure(fib.freq, "epi")
imm.tests <- cell_proportions_procedure(imm.freq, "epi")

## Epithelial compartment
tmp <- epi.freq[rownames(epi.tests$pct),colnames(epi.tests$pct)]
p3 <- matrix_barplotx(as.matrix(epi.tests$pct), tmp, group_by=epi.tests$cov$procedure, pvals=epi.tests$qvals[grepl("^resection", rownames(epi.tests$qvals)),], 
		xlab="", ylab="Percentage of sample",  colors=set.colors, legend.show=T, sig_only=F, sig_cut=0.05) + ggtitle("Epithelial cells")
pdf("Cell_composition_barplot_biopsy_vs_resection_epi6.pdf", 15, 5)
print(p3)
dev.off()
p3 <- matrix_barplotx(as.matrix(epi.tests$pct), tmp, group_by=epi.tests$cov$procedure, pvals=epi.tests$qvals[grepl("^resection", rownames(epi.tests$qvals)),], 
		xlab="", ylab="Percentage of sample",  colors=set.colors, legend.show=T, sig_only=T, sig_cut=0.05) + ggtitle("Epithelial cells")

## Stromal compartment
tmp <- fib.freq[rownames(fib.tests$pct),colnames(fib.tests$pct)]
p2 <- matrix_barplotx(as.matrix(fib.tests$pct), tmp, group_by=fib.tests$cov$procedure, pvals=fib.tests$qvals[grepl("^resection", rownames(fib.tests$qvals)),], 
		xlab="", ylab="Percentage of sample", colors=set.colors, legend.show=F, sig_only=F, sig_cut=0.05) + ggtitle("Stromal cells")
pdf("Cell_composition_barplot_biopsy_vs_resection_str6.pdf", 11, 4.4)
print(p2)
dev.off()
p2 <- matrix_barplotx(as.matrix(fib.tests$pct), tmp, group_by=fib.tests$cov$procedure, pvals=fib.tests$qvals[grepl("^resection", rownames(fib.tests$qvals)),], 
		xlab="", ylab="Percentage of sample", colors=set.colors, legend.show=F, sig_only=T, sig_cut=0.05) + ggtitle("Stromal cells")

## Immune compartment
tmp <- imm.freq[rownames(imm.tests$pct),colnames(imm.tests$pct)]
p1 <- matrix_barplotx(as.matrix(imm.tests$pct), tmp, group_by=imm.tests$cov$procedure, pvals=imm.tests$qvals[grepl("^resection", rownames(imm.tests$qvals)),], 
		xlab="", ylab="Percentage of sample", colors=set.colors, legend.show=F, sig_only=F, sig_cut=0.05) + ggtitle("Immune cells")
pdf("Cell_composition_barplot_biopsy_vs_resection_imm6.pdf", 25, 4.4)
print(p1)
dev.off()
p1 <- matrix_barplotx(as.matrix(imm.tests$pct), tmp, group_by=imm.tests$cov$procedure, pvals=imm.tests$qvals[grepl("^resection", rownames(imm.tests$qvals)),], 
		xlab="", ylab="Percentage of sample", colors=set.colors, legend.show=F, sig_only=T, sig_cut=0.05) + ggtitle("Immune cells")

pdf("Cell_composition_barplot_biopsy_vs_resection_all_sig6.pdf", 13, 4)
plot_grid(p1, p2, p3, rel_widths = c(3,6,4), nrow=1)
dev.off()
