############################################################
# Figure 3A-B  IAF score calculation and correlated cell types
############################################################

library(Seurat)
library(ggplot2)
library(cowplot)
library("RColorBrewer")
library("viridis")
library(magrittr)
library(plyr)
library(dplyr)
library(ggrepel)
library(reshape2)

seur <- readRDS(file="./data/scRNA_str.rds")
dim(seur) # 33551 15372

# Exclude inflamed cells
Idents(seur) <- seur@meta.data$status
seur <- subset(seur, idents=c("F","N","H"))
dim(seur) #  33551 14618

# Exclude epi samples
Idents(seur) <- seur@meta.data$fraction
seur <- subset(seur, idents="imu")
# subset to IAFs
Idents(seur) <- seur@meta.data$annotation2
seur <- subset(seur, idents=c("CXCL1-hi fibroblasts"))


de.cmb <- readRDS("de.cmb.rds")
de.cmb <- de.cmb[de.cmb$X=="Stromal.CXCL1-hi fibroblasts" & de.cmb$contrast=="TypeN",]
de.cmb.core <- de.cmb[de.cmb$FDR.D<0.05 & de.cmb$coefD>0,]
de.cmb.core <- de.cmb.core[order(de.cmb.core$FDR.D,decreasing=F),]
# take top 100 DEGs
top_genes <- de.cmb.core$primerid[1:100]

meta <- seur@meta.data
cpm_val <- as.matrix(GetAssayData(seur, assay = "RNA", slot = "data")[top_genes,])
dim(cpm_val) #100 1301

cpm_val <- t(cpm_val[,rownames(meta)]) # subset to paired samples
cpm_val_st <- scale(cpm_val)

fib_score <- data.frame(score=rowMeans(cpm_val_st)) %>%
	group_by(meta$sampleid, meta$patient,meta$status) %>%
	dplyr::summarise_all(mean) %>%
	as.data.frame

colnames(fib_score) <- c('sampleid', 'patient', 'status', 'score')
saveRDS(fib_score, "fib_score.rds")

sc <- fib_score$score[fib_score$status == "F"]

# add the samples dont have paired F sample
miss_s <- setdiff(unique(fib_score$patient), unique(fib_score$patient[fib_score$status == "F"][order(sc)]))

sc_ord <- c(unique(fib_score$patient[fib_score$status == "F"][order(sc)]), miss_s) 
sc_ord <- c(105387, 143371, 161862, 105446, 115284, 112584, 182097, 110201, 117351, 144119, 103292, 178961, 110645, 
			100629, 110748, 100397, 197414, 105569, 127643, 113916, 123696, 110763, 111891, 181961)

fib_score$patient <- factor(fib_score$patient, levels=sc_ord)
fib_score$patient_status <- paste(fib_score$patient, fib_score$status, sep="_")

p <- ggplot(fib_score, aes(x=factor(patient), y=score, fill=status)) + theme_bw() + 
	theme(axis.text.x=element_text(angle=45, hjust=1))+ #ggtitle("Stricture expression signature score") +
	scale_fill_manual(values=c(N="#628cd9", F="#d9628c", H="#95d962")) +
	geom_point(shape=21, size=4, alpha = 0.6) + ylab("IAF score") + xlab("DonorID")

pdf("Stricture_score1_with_heal.pdf", 8.75, 3.73)
print(p)
dev.off()


#########################################
## load cell type composition per patient

df <- read.csv("./data/updated_SS_meta.csv")
df$patient_status <- paste(df$patient, df$status, sep="_")

# keep samples that are fibrotic and imm
keep_samples <- fib_score$sampleid[fib_score$status=="F"]
df <- df[df$sampleid%in%keep_samples,]
df <- df[df$fraction=="imu",]
dim(df) # 80597    42

propmat1 <- dcast(df, annotation2 ~ sampleid)
propmat <- matrix(as.numeric(unlist(propmat1[,-1])), nrow=nrow(propmat1)) # convert to numeric
rownames(propmat) <- propmat1$annotation2
colnames(propmat) <- colnames(propmat1)[-1]

propmat <- sweep(propmat, 2, colSums(propmat), FUN="/")
fib_score.s <- fib_score[fib_score$status=="F",]
propmat_fib <- propmat[, match(fib_score.s$sampleid, colnames(propmat))]
dim(propmat_fib) # 72 21


cor.mat <- matrix(NA, nrow(propmat_fib), 2)
rownames(cor.mat) <- rownames(propmat_fib)
colnames(cor.mat) <- c("estimate", "p.val")
for(i in 1:nrow(propmat_fib)){
	res <- cor.test(fib_score.s$score, propmat_fib[i,], method="pearson")
	cor.mat[i,] <- c(res$estimate, res$p.value)
}
p.adj <- p.adjust(cor.mat[,2], method="fdr")
cor.mat <- cbind(cor.mat, p.adj)


cy <- "GZMB-hi B cells"
to_show <- bquote(atop("Pearson rho="~.(round(cor.mat[cy,1], digits = 2)), "p"["adj"]<0.05))
p1 <- ggplot(data=data.frame(x=fib_score.s$score, y=propmat_fib[rownames(propmat_fib)==cy,]), aes(x=x, y=y)) +
		annotate("text", x=0, y=0.008, label=to_show, size = unit(3, "pt")) +
		xlab("IAF score") + ylab(cy) + geom_point(color="black", size=2) + geom_smooth(method='lm', formula= y~x) + theme_bw()

cy <- "Cycling T cells"
to_show <- bquote(atop("Pearson rho="~.(round(cor.mat[cy,1], digits = 2)), "p"["adj"]<0.1))
p2 <- ggplot(data=data.frame(x=fib_score.s$score, y=propmat_fib[rownames(propmat_fib)==cy,]), aes(x=x, y=y)) +
		annotate("text", x=-0.25, y=0.01, label=to_show, size = unit(3, "pt")) +
		xlab("IAF score") + ylab(cy) + geom_point(color="black", size=2) + geom_smooth(method='lm', formula= y~x) + theme_bw()

cy <- "Follicular helper T cells"
to_show <- bquote(atop("Pearson rho="~.(round(cor.mat[cy,1], digits = 2)), "p"["adj"]<0.1))
p3 <- ggplot(data=data.frame(x=fib_score.s$score, y=propmat_fib[rownames(propmat_fib)==cy,]), aes(x=x, y=y)) +
		annotate("text", x=0, y=0.05, label=to_show, size = unit(3, "pt")) +
		xlab("IAF score") + ylab(cy) + geom_point(color="black", size=2) + geom_smooth(method='lm', formula= y~x) + theme_bw()

pdf("Stricture_score_hi_correlated_celltypes_with_heal.pdf", 9.2, 2.8)
plot_grid(p1, p2, p3, nrow=1)
dev.off()


############################################################
# Figure 3C  Barplot of genes that are correlated with IAF 
# score 
############################################################

library(Seurat)
library(ggplot2)
library(cowplot)
library("RColorBrewer")
library("viridis")
library(magrittr)
library(plyr)
library(dplyr)
library(ggrepel)
library(reshape2)
library(psych)

fib_score <- readRDS("./data/fib_score.rds")
rownames(fib_score) <- fib_score$sampleid

new_names <- read.delim("./data/Anno_names_match.txt")
dim(new_names) # 77 2

# readin different compartments
seur1 <- readRDS(file="./data/scRNA_str.rds") # stromal
seur2 <- readRDS(file="./data/scRNA_imm.rds") # immune
seur3 <- readRDS(file="./data/scRNA_epi.rds") # epithelial

# combine all comparts togetehr
seur <- merge(seur1, seur2)
seur <- merge(seur, seur3)

rm(seur1, seur2, seur3);gc()

# update the annotation2
tmp <- seur@meta.data
tmp$order <- rownames(tmp)
tmp2 <- merge(tmp, new_names, by.x="annotation2", by.y="OldName", all.x=T)
tmp2 <- tmp2[match(tmp$order, tmp2$order),]
seur$annotation2 <- tmp2$NewName


seur$annotation2_sampleid <- paste(seur$annotation2, seur$sampleid, sep=".")
Idents(seur) <- seur$annotation2_sampleid
exp <- AverageExpression(seur)$RNA
cy <- unique(seur$annotation2)
cy <- cy[grep("^Contaminant", cy, invert=T)] # exlude "Contaminant"
cy

rm(seur);gc()

##############################
# calculate the correlation, this will take some time

genescore_r <- c()
genescore_p <- c()
genescore_q <- c()
channel_num <- c()

for (i in seq_along(cy)) {
	ct_name <- cy[i]
	exp_ct <- exp[,grepl(paste("^", ct_name, sep=''), colnames(exp)), drop=F]
	dim(exp_ct)
	colnames(exp_ct) <- gsub("-", "_", gsub("^.*\\.(.*)$", "\\1", colnames(exp_ct)), fixed=T)
	exp_ct <- t(exp_ct)

	common_samples <- intersect(fib_score$sampleid, rownames(exp_ct))
	exp_ct <- exp_ct[common_samples,]
	fib_score_ct <- fib_score[common_samples,]
	
	cat(sprintf("%s %d\n", ct_name, nrow(fib_score_ct)))
	flush.console()
	
	channel_num <- c(channel_num, nrow(fib_score_ct))
	
	if(nrow(fib_score_ct)>=5){ # At least 5 channels
		exp_ct <- as.matrix(exp_ct)
		cor_ct <- corr.test(fib_score_ct[, c("score")], exp_ct)

		cor_ct$q <- matrix(p.adjust(cor_ct$p, method="fdr"), ncol=ncol(cor_ct$p))
		colnames(cor_ct$q) <- colnames(cor_ct$p)
		rownames(cor_ct$q) <- c(ct_name)
		rownames(cor_ct$r) <- c(ct_name)
		rownames(cor_ct$p) <- c(ct_name)
		
		genescore_r <- rbind(genescore_r, cor_ct$r)
		genescore_p <- rbind(genescore_p, cor_ct$p)
		genescore_q <- rbind(genescore_q, cor_ct$q)
		
	}
	gc()
}

names(channel_num) <- cy
saveRDS(list(r=genescore_r, p=genescore_p, q=genescore_q, chn=channel_num), "all_cor.scoresx.rds")

#############################################################
#  process correlation results and plot results as barplot
#############################################################

library(Seurat)
library(ggplot2)
library(cowplot)
library("RColorBrewer")
library("viridis")
library(pheatmap)
library(ggrepel)

fib_score <- readRDS("./data/fib_score.rds")
rownames(fib_score) <- fib_score$sampleid

allcor <- readRDS("./data/all_cor.scoresx.rds")
rnk <- rowSums(allcor$q < 0.05, na.rm=T)
df <- data.frame(cy=names(rnk), rnk=rnk)
dim(df) #60  2

# gather total channels per cell type
df2 <- data.frame(cy=names(allcor$chn), Channels=allcor$chn)
dim(df2) #68  2

# merge with df
df <- merge(df, df2, by="cy", all.x=T)
dim(df)  #60  3

meta <- read.csv("./data/updated_SS_meta.csv")
dim(meta) #352697     42
meta <- meta[grep("^Contaminant", meta$annotation2, invert=T),] # exlude "Contaminant"
dim(meta) # 347017     42
meta <- meta[,c("annotation2", "compart")]
meta <- unique(meta)
dim(meta) #68  2
# add compartment info
df <- merge(df, meta, by.x="cy", by.y="annotation2", all.x=T)
dim(df)  #60  4

# save all results
write.csv(df, "fibrotic_cor_scores_stats.csv")

# only take celltypes that have at least 30 channels and have at least 1 rank
df <- df[df$Channels>=30 & df$rnk>0,]
dim(df) #19  4
df <- df[df$cy!="Inflammatory fibroblasts",]
df <- df[order(df$rnk),]
df$cy <- factor(df$cy, levels=rev(df$cy))

pdf("fibrotic_cor_scores_sum.pdf", 6, 4.8)
ggplot(df, aes(y = rnk, x = cy)) + #guides(fill="none") +
  geom_bar(aes(fill = compart), stat = "identity", width = 0.5, color = "black") + theme_cowplot() +
  scale_fill_manual(values=c("Immune cells"="#0f8185", "Stromal cells"="#e86705", "Epithelial cells"="#bd409b")) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=60, vjust=1, hjust=1)) +
  ylab("Number of correlated genes") + xlab(NULL)
dev.off()

############################################################
# Figure 3F  Volcano plot of genes that are significantly
# correlatedwith IAF score in arterial endothelium 
############################################################

cy_k <- "Arterial endothelium"
keys <- allcor$q[cy_k,]
keys <- keys[!is.na(keys)]
keys_name <- names(keys)
cres <- cbind(allcor$r[cy_k, keys_name], allcor$p[cy_k, keys_name], allcor$q[cy_k, keys_name])
colnames(cres) <- c("cor","p.val","fdr")
dim(cres) # 21233     3

cres <- as.data.frame(cres)
cres$direction <- FALSE
cres$direction[cres$fdr<0.05 & cres$cor > 0] <- "POS"
cres$direction[cres$fdr<0.05 & cres$cor < 0] <- "NEG"
cres$labelx <- ""
to_show <- c("SELP", "SELL", "LAMP3", "IER5", "LIMCH1", "COL15A1", "LYVE1", "PRCP", "COL4A1", "COL4A2", "COL18A1", "COL6A2", "P4HA3", "MMP9", "HSPG2")
length(to_show) # 15
cres[to_show,]$labelx <- to_show

options(ggrepel.max.overlaps = Inf)
pdf("fibrotic_cor_scores_volcano_arterial.v6.pdf", 9.3, 7.5)
ggplot(data=cres) + geom_point(aes(x=cor, y=-log10(p.val), colour=direction), size=2) +
		scale_color_manual(values=c("POS"="#ff4040", "NEG"="#4287f5", "FALSE"="grey50")) +
		theme_bw() + #theme(legend.position="none") +
		xlab("Correlation Value") + ylab("-log10(P Value)") + 
		geom_text_repel(min.segment.length = 0, seed = 42, box.padding = 0.5, aes(x=cor, y=-log10(p.val), label=labelx))
dev.off()

############################################################
# Figure 3D  Heatmap of DEGs that are also significantly
# correlated with IAF score and their scatter plots
############################################################

library(Seurat)
library(ggplot2)
library(cowplot)
library("RColorBrewer")
library("viridis")
library(pheatmap)
library(ggrepel)

fib_score <- readRDS("./data/fib_score.rds")
rownames(fib_score) <- fib_score$sampleid

allcor <- readRDS("./data/all_cor.scoresx.rds")
de.cmb <- readRDS("./data/fib_test_de.cmb.rds")
de.cmb.o <- readRDS("de.cmb.rds")

cy_k <- "Plasma cells-IgG"
keys <- allcor$q[cy_k,]
keys <- keys[!is.na(keys)]
keys_name <- names(keys[keys<0.05])

length(de.cmb[de.cmb$cy==cy_k& de.cmb$FDR.D<0.05,]$primerid) #107
length(de.cmb[de.cmb.o$cy==cy_k& de.cmb.o$FDR.D<0.05 & de.cmb.o$contrast=="TypeN",]$primerid) #25

#################### Heatmap
convt2num.mat <- function(def, colname, type=TRUE){
	tst <- reshape2::dcast(data=def, formula = primerid~cy, fun.aggregate=NULL, value.var=colname) 
	tst.mat <- matrix(as.numeric(unlist(tst[,-1])), nrow=nrow(tst)) # convert to numeric
	rownames(tst.mat) <- tst$primerid
	colnames(tst.mat) <- colnames(tst)[-1]
	
	if(type){
		tst.mat[is.na(tst.mat)] <- 0 #change NAs to 0
	}else{ # fdr or p.value
		tst.mat[is.na(tst.mat)] <- 1 #change NAs to 1
	}

	return(tst.mat)
}

## 
most_deg.f <- de.cmb[de.cmb$cy==cy_k& de.cmb$FDR.D<0.05,]
keys1 <- intersect(most_deg.f$primerid, keys_name)
length(keys1)
most_deg.f <- de.cmb.o[de.cmb.o$cy==cy_k& de.cmb.o$FDR.D<0.05 & de.cmb.o$contrast=="TypeN",]
keys2 <- intersect(most_deg.f$primerid, keys_name)
length(keys2)
length(intersect(keys1,keys2))
length(union(keys1,keys2))
comm <- union(keys1,keys2)

most_deg.f <- de.cmb[de.cmb$cy==cy_k,]
most_deg.f <- most_deg.f[most_deg.f$primerid%in%comm, ]
num.mat1 <- convt2num.mat(most_deg.f, "coefD")
fdr.mat1 <- convt2num.mat(most_deg.f, "FDR.D", type=FALSE)
num.mat1 <- num.mat1[comm,]
fdr.mat1 <- fdr.mat1[comm,]

most_deg.f <- de.cmb.o[de.cmb.o$cy==cy_k & de.cmb.o$contrast=="TypeN",]
most_deg.f <- most_deg.f[most_deg.f$primerid%in%comm, ]
num.mat2 <- convt2num.mat(most_deg.f, "coefD")
fdr.mat2 <- convt2num.mat(most_deg.f, "FDR.D", type=FALSE)
num.mat2 <- num.mat2[comm,]
fdr.mat2 <- fdr.mat2[comm,]

num.mat <- cbind(num.mat1, num.mat2)
colnames(num.mat) <- c("IAF","All")
fdr.mat <- cbind(fdr.mat1, fdr.mat2)
pdf("fibrotic_cor_scores_heatmap_plasma_igg.pdf", 4.5, 4.8)
pheatmap(t(num.mat), color=colorRampPalette(c("#0068bd", "white", "#bd0000"))(100), cluster_cols=T, cluster_rows=F,
		fontsize = 12, show_rownames=T, show_colnames=T, angle_col = "315", breaks=seq(-2,2,length=101), cellwidth=20, cellheight=20,
		display_numbers = matrix(ifelse(t(fdr.mat)<0.05, "*", ""), nrow(t(fdr.mat))), number_color="black", border_color = "grey")
dev.off()

#############################################
## plot scatter

seur <- readRDS(file="./data/ss_imm.rds")
dim(seur) # 33551 156147

# Exclude UC patient 174879
Idents(seur) <- seur@meta.data$patient
seur <- subset(seur, idents=174879, invert = TRUE)
dim(seur) #33551 150759

Idents(seur) <- seur@meta.data$clean
seur <- subset(seur, idents="no", invert = TRUE)
dim(seur) #33551 150754

cy_k #"Plasma cells-IgG"
Idents(seur) <- seur$annotation2
seur <- subset(seur, idents=cy_k)
dim(seur)
seur$annotation2_sampleid <- paste(seur$annotation2, seur$sampleid, sep=".")
Idents(seur) <- seur$annotation2_sampleid
exp_ct <- AverageExpression(seur)$RNA
colnames(exp_ct) <- gsub("-", "_", gsub("^.*\\.(.*)$", "\\1", colnames(exp_ct)), fixed=T)
exp_ct <- t(exp_ct)

common_samples <- intersect(fib_score$sampleid, rownames(exp_ct))
exp_ct <- exp_ct[common_samples,]
fib_score_ct <- fib_score[common_samples,]
exp_ct <- as.matrix(exp_ct)

comm #"GAS6"      "MAP4K3-DT" "PTP4A3"    "TNFAIP2"  
allcor$r[cy_k,comm]

y_pos <- c(1, 0.3, 0.95, 0.9)
plotlist <- list()
for(i in 1:length(comm)){ # top pairs
	#to_show <- sprintf("Pearson\n rho=%s",round(allcor$r[cy_k,comm[i]], digits = 2))
	to_show <- bquote(atop("Pearson rho="~.(round(allcor$r[cy_k,comm[i]], digits = 2)), "p"["adj"]<0.05))
	p1 <- ggplot(data=data.frame(x=fib_score_ct$score, y=exp_ct[,comm[i]]), aes(x=x, y=y)) +
		geom_point() + theme_bw() +
		geom_smooth(method="lm", formula=y~x) +
		xlab("IAF score") + ylab(comm[i]) +
		annotate("text", x=0, y=y_pos[i], label=to_show, size = unit(3, "pt"))
		plotlist <- c(plotlist, list(p1))
}
length(plotlist)

pdf("fibrotic_cor_scores_scatters_plasma_igg.v6.1.pdf", 5.6, 5.6)
plot_grid(plotlist=plotlist, nrow=2, ncol=2)
dev.off()
