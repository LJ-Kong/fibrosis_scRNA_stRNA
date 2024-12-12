############################################################
# Figure 2A   Compositional tests
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
meta <- read.csv("updated_SS_meta.csv")
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

epi.tests <- cell_proportions_tests(epi.freq, "imu")
fib.tests <- cell_proportions_tests(fib.freq, "epi")
imm.tests <- cell_proportions_tests(imm.freq, "epi")


## Epithelial compartment
tmp <- epi.freq[rownames(epi.tests$pct),colnames(epi.tests$pct)]
p3 <- matrix_barplot(as.matrix(epi.tests$pct), tmp, group_by=epi.tests$cov$condition, pvals=epi.tests$qvals[3:4,], 
		xlab="", ylab="Percentage of sample",  colors=set.colors, legend.show=T, sig_only=F, sig_cut=0.05) + ggtitle("Epithelial cells")
pdf("Cell_composition_barplot_epi6.pdf", 15, 5)
print(p3)
dev.off()
# only significant
p3 <- matrix_barplot(as.matrix(epi.tests$pct), tmp, group_by=epi.tests$cov$condition, pvals=epi.tests$qvals[3:4,], 
		xlab="", ylab="Percentage of sample", colors=set.colors, legend.show=T, sig_only=T, sig_cut=0.05) + ggtitle("Epithelial cells")
pdf("Cell_composition_barplot_epi_only_sig6.pdf", 4.8, 4)
print(p3)
dev.off()


## Stromal compartment
tmp <- fib.freq[rownames(fib.tests$pct),colnames(fib.tests$pct)]
p2 <- matrix_barplot(as.matrix(fib.tests$pct), tmp, group_by=fib.tests$cov$condition, pvals=fib.tests$qvals[3:4,], 
		xlab="", ylab="Percentage of sample", colors=set.colors, legend.show=F, sig_only=F, sig_cut=0.05) + ggtitle("Stromal cells")
pdf("Cell_composition_barplot_str6.pdf", 11, 4.4)
print(p2)
dev.off()
# only significant
p2 <- matrix_barplot(as.matrix(fib.tests$pct), tmp, group_by=fib.tests$cov$condition, pvals=fib.tests$qvals[3:4,], 
		xlab="", ylab="Percentage of sample", colors=set.colors, legend.show=T, sig_only=T, sig_cut=0.05) + ggtitle("Stromal cells")
pdf("Cell_composition_barplot_str_only_sig6.pdf", 5.9, 4)
print(p2)
dev.off()


## Immune compartment
tmp <- imm.freq[rownames(imm.tests$pct),colnames(imm.tests$pct)]
p1 <- matrix_barplot(as.matrix(imm.tests$pct), tmp, group_by=imm.tests$cov$condition, pvals=imm.tests$qvals[3:4,], 
		xlab="", ylab="Percentage of sample", colors=set.colors, legend.show=F, sig_only=F, sig_cut=0.05) + ggtitle("Immune cells")
pdf("Cell_composition_barplot_imm6.pdf", 25, 4.4)
print(p1)
dev.off()
# only significant
p1 <- matrix_barplot(as.matrix(imm.tests$pct), tmp, group_by=imm.tests$cov$condition, pvals=imm.tests$qvals[3:4,], 
		xlab="", ylab="Percentage of sample", colors=set.colors, legend.show=F, sig_only=T, sig_cut=0.05) + ggtitle("Immune cells")
pdf("Cell_composition_barplot_imm_only_sig6.pdf", 8, 4)
print(p1)
dev.off()


# Count immune cells that are from Epithelial samples
immx.freq = as.matrix(as.data.frame.matrix(table(meta$sampleid[meta$compart=="Immune cells"& meta$fraction=="epi"], meta$anno22[meta$compart=="Immune cells"& meta$fraction=="epi"])))
immx.tests <- cell_proportions_tests(immx.freq, "imu")

tmp <- immx.freq[rownames(immx.tests$pct),colnames(immx.tests$pct)]
p1 <- matrix_barplot(as.matrix(immx.tests$pct), tmp, group_by=immx.tests$cov$condition, pvals=immx.tests$qvals[3:4,], 
		xlab="", ylab="Percentage of sample", colors=set.colors, legend.show=F, sig_only=F, sig_cut=0.05) + ggtitle("Epithelial fraction immune cells")
pdf("Cell_composition_barplot_imm6_intraepi.pdf", 25, 4.4)
print(p1)
dev.off()
# only significant
p1 <- matrix_barplot(as.matrix(immx.tests$pct), tmp, group_by=immx.tests$cov$condition, pvals=immx.tests$qvals[3:4,], 
		xlab="", ylab="Percentage of sample", colors=set.colors, legend.show=F, sig_only=T, sig_cut=0.05) + ggtitle("Epithelial fraction immune cells")
pdf("Cell_composition_barplot_imm6_intraepi_only_sig.pdf", 2.8, 4)
print(p1)
dev.off()


############################################################
# Figure 2B   Paired compositional tests
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
library(brms)

source("dirichlet_with_random.r")
source("dirichlet_functions.r")

# Load metadata
meta <- read.csv("updated_SS_meta.csv")

# get pair info
df2 <- read.csv("combine_samples_by_fraction_manual3.csv")
dim(df2) # 102   9

# only keep the paired epi/imu samples
meta <- meta[meta$sampleid %in% df2$sampleid[df2$paired_epi=="Y" | df2$paired_imu=="Y"],]
dim(meta) # 219865     41

meta <- meta[grep("^Contaminant", meta$annotation2, invert=T),] # exlude "Contaminant"
dim(meta) #  216191     41

# one imu sample has no stromal cells at all
table(meta$compart[meta$sampleid=="105446_N14_imu"])


###########################
# create a new column change anno2 from factor to character
meta$anno22 <- as.character(meta$annotation2)

epi.freq = as.matrix(as.data.frame.matrix(table(meta$sampleid[meta$compart=="Epithelial cells"], meta$anno22[meta$compart=="Epithelial cells"])))
fib.freq = as.matrix(as.data.frame.matrix(table(meta$sampleid[meta$compart=="Stromal cells"], meta$anno22[meta$compart=="Stromal cells"])))
fib.freq = fib.freq[!rownames(fib.freq)%in%c("105446_N14_imu","105446_F14_imu"),] # exclude the pair of samples which one of them doesnt have stromal cells
imm.freq = as.matrix(as.data.frame.matrix(table(meta$sampleid[meta$compart=="Immune cells"], meta$anno22[meta$compart=="Immune cells"])))

# epi
epi.fit <- cell_proportions_tests_random(epi.freq, "imu", iter=6000)
epi.tests <- cell_proportions_tests(epi.freq, "imu")
tmp1 <- hypothesis(epi.fit, paste(epi.fit$family$dpars[-length(epi.fit$family$dpars)], "_conditionN=0", sep=""))


pass_pvals <- as.matrix(t(tmp1$hypothesis$Post.Prob))
colnames(pass_pvals) <- sort(colnames(epi.tests$pct))
rownames(pass_pvals) <- "N"
p3 <- matrix_barplot_random(as.matrix(epi.tests$pct), epi.freq[rownames(epi.tests$pct),colnames(epi.tests$pct)], group_by=epi.tests$cov$condition, 
						pvals=pass_pvals, xlab="", ylab="",  colors=set.colors, legend.show=T, sig_only=F, sig_cut=0.05) + ggtitle("Epithelial cells")
pdf("Cell_composition_barplot_epi_bayes_paired2.pdf", 11, 4)
print(p3)
dev.off()
# only significant
p3 <- matrix_barplot_random(as.matrix(epi.tests$pct), epi.freq[rownames(epi.tests$pct),colnames(epi.tests$pct)], group_by=epi.tests$cov$condition, 
						pvals=pass_pvals, xlab="", ylab="",  colors=set.colors, legend.show=T, sig_only=T, sig_cut=0.05) + ggtitle("Epithelial cells")
pdf("Cell_composition_barplot_epi_bayes_paired_sig2.pdf", 4.86, 4)
print(p3)
dev.off()


# str
fib.fit <- cell_proportions_tests_random(fib.freq, "epi", iter=6000)
fib.tests <- cell_proportions_tests(fib.freq, "epi")
tmp2 <- hypothesis(fib.fit, paste(fib.fit$family$dpars[-length(fib.fit$family$dpars)], "_conditionN=0", sep=""))

pass_pvals <- as.matrix(t(tmp2$hypothesis$Post.Prob))
colnames(pass_pvals) <- sort(colnames(fib.tests$pct))
rownames(pass_pvals) <- "N"
p2 <- matrix_barplot_random(as.matrix(fib.tests$pct), fib.freq[rownames(fib.tests$pct),colnames(fib.tests$pct)], group_by=fib.tests$cov$condition, 
						pvals=pass_pvals, xlab="", ylab="",  colors=set.colors, legend.show=F, sig_only=F, sig_cut=0.05) + ggtitle("Stromal cells")
pdf("Cell_composition_barplot_str_bayes_paired2.pdf", 9, 4)
print(p2)
dev.off()
# only significant
p2 <- matrix_barplot_random(as.matrix(fib.tests$pct), fib.freq[rownames(fib.tests$pct),colnames(fib.tests$pct)], group_by=fib.tests$cov$condition, 
						pvals=pass_pvals, xlab="", ylab="",  colors=set.colors, legend.show=T, sig_only=T, sig_cut=0.05) + ggtitle("Stromal cells")
pdf("Cell_composition_barplot_str_bayes_paired_sig2.pdf", 4.5, 4)
print(p2)
dev.off()


imm.fit <- cell_proportions_tests_random(imm.freq, "epi", iter=6000)
imm.tests <- cell_proportions_tests(imm.freq, "epi")
tmp3 <- hypothesis(imm.fit, paste(imm.fit$family$dpars[-length(imm.fit$family$dpars)], "_conditionN=0", sep=""))

pass_pvals <- as.matrix(t(tmp3$hypothesis$Post.Prob))
colnames(pass_pvals) <- sort(colnames(imm.tests$pct))
rownames(pass_pvals) <- "N"
p1 <- matrix_barplot_random(as.matrix(imm.tests$pct), imm.freq[rownames(imm.tests$pct),colnames(imm.tests$pct)], group_by=imm.tests$cov$condition, 
						pvals=pass_pvals, xlab="", ylab="",  colors=set.colors, legend.show=F, sig_only=F, sig_cut=0.05) + ggtitle("Immune cells")
pdf("Cell_composition_barplot_imm_bayes_paired2.pdf", 21, 4)
print(p1)
dev.off()
# only significant
p1 <- matrix_barplot_random(as.matrix(imm.tests$pct), imm.freq[rownames(imm.tests$pct),colnames(imm.tests$pct)], group_by=imm.tests$cov$condition, 
						pvals=pass_pvals, xlab="", ylab="",  colors=set.colors, legend.show=T, sig_only=T, sig_cut=0.05) + ggtitle("Immune cells")
pdf("Cell_composition_barplot_imm_bayes_paired_sig2.pdf", 5.3, 4)
print(p1)
dev.off()


meta$tissue[meta$tissue=="ascending colon"] <- "colon"
immx.freq = as.matrix(as.data.frame.matrix(table(meta$sampleid[meta$compart=="Immune cells"& meta$fraction=="epi"], meta$anno22[meta$compart=="Immune cells"& meta$fraction=="epi"])))
immx.fit <- cell_proportions_tests_random(immx.freq, "imu", iter=6000)
immx.tests <- cell_proportions_tests(immx.freq, "imu")
tmp3 <- hypothesis(immx.fit, paste(immx.fit$family$dpars[-length(immx.fit$family$dpars)], "_conditionN=0", sep=""))

pass_pvals <- as.matrix(t(tmp3$hypothesis$Post.Prob))
colnames(pass_pvals) <- sort(colnames(immx.tests$pct))
rownames(pass_pvals) <- "N"
p1 <- matrix_barplot_random(as.matrix(immx.tests$pct), immx.freq[rownames(immx.tests$pct),colnames(immx.tests$pct)], group_by=immx.tests$cov$condition, 
						pvals=pass_pvals, xlab="", ylab="",  colors=set.colors, legend.show=F, sig_only=F, sig_cut=0.05) + ggtitle("Intraepithelial cells")
pdf("Cell_composition_barplot_imm_bayes_paired_intraepi2.pdf", 21, 4)
print(p1)
dev.off()
# only significant
p1 <- matrix_barplot_random(as.matrix(immx.tests$pct), immx.freq[rownames(immx.tests$pct),colnames(immx.tests$pct)], group_by=immx.tests$cov$condition, 
						pvals=pass_pvals, xlab="", ylab="",  colors=set.colors, legend.show=T, sig_only=T, sig_cut=0.05) + ggtitle("Intraepithelial cells")
pdf("Cell_composition_barplot_imm_bayes_paired_intraepi_sig2.pdf", 5.3, 4)
print(p1)
dev.off()





############################################################
# Figure 2C  Barplots of differentially expressed (DE) genes
############################################################
# Please run DE_run.r to get the DE results


# readin updated annotation for cell types
new_names <- read.delim("Anno_names_match.txt")
# load DE result 
de.cmb <- readRDS("de.cmb.rds")
# subset to significant ones
de.cmb.sub <- de.cmb[de.cmb$FDR.D<0.05,] 

# Total DEGs in F vs. H
dim(de.cmb[de.cmb$FDR.D<0.05 & de.cmb$contrast=="TypeH",])
# Total DEGs in F vs. N
dim(de.cmb[de.cmb$FDR.D<0.05 & de.cmb$contrast=="TypeN",])

df <- de.cmb.sub
df1 <- df[df$contrast=="TypeN", ]
dfgg1c <-  dplyr::count(df1, X, direct.D,.drop=FALSE)
dfgg1c$direct.D <- factor(dfgg1c$direct.D, levels=c(FALSE, TRUE), labels=c("UP", "DOWN"))
dfgg1c$loc <- gsub("^(\\w+)\\..*$", "\\1", dfgg1c$X)
dfgg1c$cy <- gsub("^\\w+\\.", "", dfgg1c$X)
dfgg1cx <- merge(dfgg1c, new_names, by.x="cy", by.y="OldName", all.x=T) #update names
dfgg1cx <- dfgg1cx[grep("^Contaminant", dfgg1cx$NewName, invert=T),] # exlude "Contaminant"


df2 <- df[df$contrast=="TypeH", ]
dfgg2c <-  dplyr::count(df2, X, direct.D,.drop=FALSE)
dfgg2c$direct.D <- factor(dfgg2c$direct.D, levels=c(FALSE, TRUE), labels=c("UP", "DOWN"))
dfgg2c$loc <- gsub("^(\\w+)\\..*$", "\\1", dfgg2c$X)
dfgg2c$cy <- gsub("^\\w+\\.", "", dfgg2c$X)
dfgg2cx <- merge(dfgg2c, new_names, by.x="cy", by.y="OldName", all.x=T) #update names
dfgg2cx <- dfgg2cx[grep("^Contaminant", dfgg2cx$NewName, invert=T),] # exlude "Contaminant"


legend_title <- "Direction"

library(grid)

p1 <- ggplot(dfgg1cx)+ labs(title="Stricture vs. Non-stricture", fill=legend_title) + theme_bw() + ylim(0, 600) + facet_grid(loc ~ ., scales="free_y",space="free") +
		theme(axis.title.y = element_blank(), strip.text.y = element_text(face="bold")) + 
		geom_bar(aes(x=NewName, y=n, fill=direct.D), stat="identity", width=0.6, position = position_stack(reverse = T))+ coord_flip() + ylab("Number of DEGs")
g <- ggplot_gtable(ggplot_build(p1))
strip_r <- which(grepl('strip-r', g$layout$name))
fills <- c("#ecd4ff", "#ffbe85", "#cf9286")
k <- 1
for (i in strip_r) {
	j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
	g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
	k <- k+1
}
pdf("DEGs_sum_barplot_noni-fibr.pdf", 6, 7.5)
grid.draw(g)
dev.off()

p2 <- ggplot(dfgg2cx)+ labs(title="Stricture vs. Non-IBD control", fill=legend_title) + theme_bw() + ylim(0, 600) + facet_grid(loc ~ ., scales="free_y",space="free") +
		theme(axis.title.y = element_blank(), strip.text.y = element_text(face="bold")) +
		geom_bar(aes(x=NewName, y=n, fill=direct.D), stat="identity", width=0.6, position = position_stack(reverse = T))+ coord_flip() + ylab("Number of DEGs")
g <- ggplot_gtable(ggplot_build(p2))
strip_r <- which(grepl('strip-r', g$layout$name))
fills <- c("#ecd4ff", "#ffbe85", "#cf9286")
k <- 1
for (i in strip_r) {
	j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
	g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
	k <- k+1
}
pdf("DEGs_sum_barplot_heal-fibr.pdf", 6, 7.5)
grid.draw(g)
dev.off()



############################################################
# Figure 2D  Scatter plots of coefD
############################################################

library(ggplot2)
library(cowplot)


# load DE result 
de.cmb <- readRDS("./data/de.cmb.rds")
dim(de.cmb) #606557     21

scatter.coef <- function (cmp, xp, yp) {
	df <- de.cmb[de.cmb$loc==cmp,]
	
	df <- df[,c("primerid", "coefD", "FDR.D", "cy", "contrast")]
	df <- df[complete.cases(df),]
	
	df <- df[abs(df$coefD)<10, ] # make a cut so there is no crazy point
	var_regex = '^MT|^RP' # remove MT, and RP genes based on HUGO gene names
	df <- df[grep(var_regex, df$primerid, invert=T),]
	
	df$cmb <- paste(df$primerid, df$cy, sep=".")
	ids <- unique(df$cmb);length(ids) #67246

	dfx <- df[df$contrast=="TypeH",]
	dfy <- df[df$contrast=="TypeN",]
	dfz <- data.frame(Fibro.vs.Heal=dfx[match(ids, dfx$cmb),]$coefD, Fibro.vs.NonF=dfy[match(ids, dfy$cmb),]$coefD)
	rownames(dfz) <- ids
	dfz <- dfz[complete.cases(dfz),]
	
	cor.res <- cor.test(dfz$Fibro.vs.Heal, dfz$Fibro.vs.NonF, method = "spearman", alternative = "greater")
	cor.info <- sprintf("Spearman cor=%f\n(p.value=%e)", cor.res$estimate, cor.res$p.value)
	cat(cor.info)
	flush.console()
	mtitle <- cmp
	
	#to_show <- sprintf("Pearson rho=%s \n p<10^-307",round(cor.res$estimate, digits = 2))
	to_show <- bquote(atop("Spearman rho ="~.(round(cor.res$estimate, digits = 2)), "p"<10^-307))

	p <- ggplot(dfz, aes(x=Fibro.vs.Heal, y=Fibro.vs.NonF)) + xlab("DE coefficient\nStricture vs. Non-IBD") + ylab("DE coefficient\nNon-stricture vs. Non-IBD") +
			geom_point(color="black", size=0.55) + theme_bw(base_size = 12) + labs(title=mtitle) +
			geom_abline(intercept = 0, slope = 1, color="tan2") +
			annotate("text", x=xp, y=yp, label=to_show, size = unit(2.5, "pt"))
			theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	
	return(p)
}

sc1 <- scatter.coef("Epithelial", -4, 3.5) 
sc2 <- scatter.coef("Stromal", -2.5, 4) 
sc3 <- scatter.coef("Immune", -4.5, 7) 

pvt <- plot_grid(sc1,sc2,sc3,nrow=1)
pdf("fig_scatter.coefD.pdf", 8.6, 3.2)
pvt
dev.off()



############################################################
# Figure 2E  Pathway analysis
############################################################


library(ggplot2)
library(cowplot)
library("RColorBrewer")
library("viridis")
library(dplyr)
library(ggrepel)
library(reshape2)
library(MASS)
library(egg)
library(fgsea)
library(pheatmap)


## Download pathways from: https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp#C2
pty <- "./data/c2.cp.kegg_legacy.v2023.2.Hs.symbols.gmt"

# load DE result 
de.cmb <- readRDS("./data/de.cmb.rds")

out_col <- c("TypeH", "TypeN")
cmb.fgseaRes <- c()
for(j in 1:2){

	# go through all celltypes 
	celltys <- unique(de.cmb$cy)

	tmp.fgseaRes <- c()
	for(i in 1:length(celltys)){
		cellty <- celltys[i]
		cat(sprintf("Working on: [%d] %s\n", i, cellty))
		flush.console()
		
		# use pre.coefD!! But the direction needs to be flipped!!
		ranklist <- -1*de.cmb$coefD.pre[de.cmb$cy==cellty & de.cmb$contrast==out_col[j]] 
		names(ranklist) <- de.cmb$primerid[de.cmb$cy==cellty & de.cmb$contrast==out_col[j]]
		ranklist <- ranklist[!is.na(ranklist)]
		
		if(length(ranklist)>100){
			fgseaRes <- fgsea(pathways = pathways.hallmark, 
								  stats = sort(ranklist),
								  minSize=3,
								  maxSize=500,
								  nperm=100000)
			
			fgseaRes$cellty <- cellty
			tmp.fgseaRes <- rbind(tmp.fgseaRes, fgseaRes)
		}
	}
	
	tmp.fgseaRes$contrast <- out_col[j]
	if(j==1){
		cmb.fgseaRes <- tmp.fgseaRes
	}else{
		cmb.fgseaRes <- rbind(cmb.fgseaRes, tmp.fgseaRes)
	}
}

dim(cmb.fgseaRes) 

cpart <- de.cmb[,c("loc","cy")]
cpart <- cpart[!duplicated(cpart),]

# add compartment info
cmb.fgseaRes <- merge(cmb.fgseaRes, cpart, by.x="cellty", by.y="cy", all.x=T)

# readin updated annotation for cell types
new_names <- read.delim("Anno_names_match.txt")
cmb <- merge(cmb.fgseaRes, new_names, by.x="cellty", by.y="OldName", all.x=T) #update names
colnames(cmb)[1] <- "OldName"
colnames(cmb)[12] <- "cellty"

#saveRDS(cmb, "cmb.fgseaRes.coefD.pre.kegg_legacy.rds")
cmb <- readRDS("cmb.fgseaRes.coefD.pre.kegg_legacy.rds")
cmb <- cmb[grep("^Contaminant", cmb$cellty, invert=T),] # exlude "Contaminant"


## pwy.collect returns 2 matrice: NES and padj (at least one significant), 
## and a data.frame contains total number of significant pwy per each direction
pwy.collect <- function(fRes, is_kegg=T){
	fRes$cellty <- paste(fRes$loc, fRes$cellty, sep=".")
	tst <- dcast(data=fRes, formula = pathway~cellty, fun.aggregate=sum, value.var="NES")
	tst <- as.data.frame(tst) # class(tst) -> "data.table" "data.frame"
	annotation_col = data.frame(CellType=substr(colnames(tst)[-1], 1, 3))
	colnames(tst) <- gsub("^\\w+\\.", "", colnames(tst))
	rownames(annotation_col) <- colnames(tst)[-1]
	if(is_kegg){
		rownames(tst) <- substr(tst$pathway, 6, nchar(tst$pathway))
	}else{
		rownames(tst) <- substr(tst$pathway, 10, nchar(tst$pathway))
	}
	tst <- tst[,-1]
	

	tst2 <- dcast(data=fRes, formula = pathway~cellty, fun.aggregate=sum, value.var="padj")
	tst2 <- as.data.frame(tst2)
	colnames(tst2) <- gsub("^\\w+\\.", "", colnames(tst2))
	if(is_kegg){
		rownames(tst2) <- substr(tst2$pathway, 6, nchar(tst2$pathway))
	}else{
		rownames(tst2) <- substr(tst2$pathway, 10, nchar(tst2$pathway))
	}
	tst2 <- tst2[,-1]
	tst2[tst2==0] <- 1 # change all 0 padj to 1

	sign_paths <- dplyr::count(fRes[fRes$padj<0.05,],pathway)
	if(is_kegg){
		sign_paths$pathway <- substr(sign_paths$pathway, 6, nchar(sign_paths$pathway))
	}else{
		sign_paths$pathway <- substr(sign_paths$pathway, 10, nchar(sign_paths$pathway))
	}
	sign_paths <- sign_paths[order(sign_paths$n, decreasing=T),]
	
	sign_tst <- tst[sign_paths$pathway,]
	sign_tst <- sign(sign_tst)
	
	sign_tst2 <- tst2[sign_paths$pathway,]
	sign_tst2[sign_tst2<0.05] <- -1
	sign_tst2[sign_tst2>=0.05] <- 0
	sign_tst2[sign_tst2==-1] <- 1
	
	sign_cmb <- sign_tst*sign_tst2 # collect NES directions and padj together
	sign_cmb_stat <- data.frame(pos=apply(sign_cmb, 1, function(x) sum(x>0)), neg=apply(sign_cmb, 1, function(x) sum(x<0)))
	rownames(sign_cmb_stat) <- rownames(sign_cmb)
	
	return(list(nes.mat=tst, padj.mat=tst2, collect.stat=sign_cmb_stat, anno_col_df=annotation_col))
}

########################################################
## Epi
fgseaRes <- cmb[cmb$contrast=="TypeH"&cmb$loc=="Epithelial",] #TypeH, Epi
res1 <- pwy.collect(fgseaRes)


fgseaRes <- cmb[cmb$contrast=="TypeN"&cmb$loc=="Epithelial",] #TypeN, Epi
res2 <- pwy.collect(fgseaRes)

## Str
fgseaRes <- cmb[cmb$contrast=="TypeH"&cmb$loc=="Stromal",] #TypeH, Str
res3 <- pwy.collect(fgseaRes)

fgseaRes <- cmb[cmb$contrast=="TypeN"&cmb$loc=="Stromal",] #TypeN, Str
res4 <- pwy.collect(fgseaRes)

## Imm
fgseaRes <- cmb[cmb$contrast=="TypeH"&cmb$loc=="Immune",] #TypeH, Imm
res5 <- pwy.collect(fgseaRes)

fgseaRes <- cmb[cmb$contrast=="TypeN"&cmb$loc=="Immune",] #TypeN, Imm
res6 <- pwy.collect(fgseaRes)


###############
## Epi, barplot
df1 <- res1$collect.stat
df1$pwy <- rownames(df1)
df1$loc <- "H" 
df2 <- res2$collect.stat
df2$pwy <- rownames(df2)
df2$loc <- "N" 
df <- rbind(df1,df2)
df <- melt(df)

df$pct[df$loc=="H"] <- df$value[df$loc=="H"]/ncol(res1$nes.mat)
df$pct[df$loc=="N"] <- df$value[df$loc=="N"]/ncol(res2$nes.mat)
df$loc <- factor(df$loc, levels=c("H","N"))
df$cpart <- "Epi"
df.epi <- df


###############
## Str, barplot
df1 <- res3$collect.stat
df1$pwy <- rownames(df1)
df1$loc <- "H" 
df2 <- res4$collect.stat
df2$pwy <- rownames(df2)
df2$loc <- "N" 
df <- rbind(df1,df2)
df <- melt(df)

df$pct[df$loc=="H"] <- df$value[df$loc=="H"]/ncol(res3$nes.mat)
df$pct[df$loc=="N"] <- df$value[df$loc=="N"]/ncol(res4$nes.mat)
df$loc <- factor(df$loc, levels=c("H","N"))
df$cpart <- "Str"
df.str <- df


###############
## Imm, barplot
df1 <- res5$collect.stat
df1$pwy <- rownames(df1)
df1$loc <- "H" 
df2 <- res6$collect.stat
df2$pwy <- rownames(df2)
df2$loc <- "N" 
df <- rbind(df1,df2)
df <- melt(df)

df$pct[df$loc=="H"] <- df$value[df$loc=="H"]/ncol(res5$nes.mat)
df$pct[df$loc=="N"] <- df$value[df$loc=="N"]/ncol(res6$nes.mat)
df$loc <- factor(df$loc, levels=c("H","N"))
df$cpart <- "Imm"
df.imm <- df

#############
df.cmb <- rbind(df.epi, df.str, df.imm)

pct.mat <- dcast(data=df.cmb, formula = pwy~cpart+loc, fun.aggregate=sum, value.var="pct")
pct.mat <- pct.mat[rowSums(pct.mat[,-1]>0.1)>1,] # filter out some low percentage pwy
pct.mat.f <- dcast(data=df.cmb, formula = pwy~variable+cpart+loc, fun.aggregate=sum, value.var="pct") # add direction back
pct.mat.f <- pct.mat.f[pct.mat.f$pwy%in%pct.mat$pwy,] # keep the useful pwy
dim(pct.mat.f) #33  13

hc <- hclust(dist(pct.mat.f[,-1])) # order the pct
order_pwy <- pct.mat.f$pwy[hc$order]

df.pct.mat.f <- melt(pct.mat.f)
df.pct.mat.f$direction <- factor(substr(df.pct.mat.f$variable, 1, 3), levels=c("neg", "nil", "pos"))
df.pct.mat.f$cpart <- substr(df.pct.mat.f$variable, 5, 7)
df.pct.mat.f$loc <- substr(df.pct.mat.f$variable, 9, 10)
df.pct.mat.f$loc <- factor(df.pct.mat.f$loc, levels=c("H","N"), labels=c("Stricture vs. Non-IBD","Stricture vs. Non-stricture"))
df.pct.mat.f$pwy <- factor(df.pct.mat.f$pwy, levels=order_pwy)
df.pct.mat.f$cpart[df.pct.mat.f$cpart=="Epi"] <- "Epithelial"
df.pct.mat.f$cpart[df.pct.mat.f$cpart=="Str"] <- "Stromal"
df.pct.mat.f$cpart[df.pct.mat.f$cpart=="Imm"] <- "Immune"

# kegg_legacy
pdf("pathways_barplot_sum_kegg_legacy.pdf", 16, 5)
ggplot(data=df.pct.mat.f) + 
		geom_bar(aes(x=pwy, y=ifelse(direction=="neg",-value, value), fill=direction), stat="identity",
		width=0.6, position = position_stack(reverse = T))+ theme_bw() +
		facet_grid(~cpart+loc) + coord_flip() + 
		#theme(axis.text.y = element_text(color = "grey20", size = 5)) +
		scale_fill_manual(name="Enrichment direction", labels=c(pos="Positive", neg="Negative"), values=c(pos="tomato3", neg="turquoise4")) +
		 ylab("Fraction of cell types") + xlab("")
dev.off()


