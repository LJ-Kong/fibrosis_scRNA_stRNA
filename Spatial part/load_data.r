
# ------------------------------------------------------------
# IMPORTANT: replace this with the path to your project folder
# ------------------------------------------------------------
PROJECT_FOLDER = '/broad/smillie-data/proj/cd-spatial/github/'
if(! file.exists(PROJECT_FOLDER)){
    stop('PROJECT_FOLDER does not exist. You need to modify the path at top of load_data.r')
}

get_path = function(path, base_dir=PROJECT_FOLDER){paste0(gsub('/+$', '', base_dir), '/', path)}

# ---------
# load code
# ---------
library(akima)
library(biomaRt)
library(BiocGenerics)
library(cccd)
library(colorspace)
library(cowplot)
library(data.table)
library(DelayedArray)
library(entropy)
library(expm)
library(FNN)
library(ggbeeswarm)
library(ggdendro)
library(ggnetwork)
library(ggplot2)
library(ggrepel)
library(ggsignif)
library(gplots)
library(grDevices)
library(grid)
library(gtools)
library(harmony)
library(Hmisc)
library(hues)
library(igraph)
library(irlba)
library(lme4)
library(MASS)
library(MAST)
library(Matrix)
library(Matrix.utils)
library(methods)
library(mygene)
library(naturalsort)
library(NMF)
library(org.Hs.eg.db)
library(png)
library(proxy)
library(RColorBrewer)
library(rentrez)
library(repr)
library(rhdf5)
library(rliger)
library(rjson)
library(rsvd)
library(S4Vectors)
library(scales)
library(sva)
library(ROCR)
library(Rtsne)
library(tibble)
library(tidyr)
library(tidyverse)
library(umap)
library(viridis)
library(wordspace)
library(zoo)


source(get_path('code/batch.r'))
source(get_path('code/cluster.r'))
source(get_path('code/colors.r'))
source(get_path('code/dendro.r'))
source(get_path('code/downsample.r'))
source(get_path('code/gene_module.r'))
source(get_path('code/map_gene.r'))
source(get_path('code/markers.r'))
source(get_path('code/mm_utils.r'))
source(get_path('code/mtx.r'))
source(get_path('code/network.r'))
source(get_path('code/pca.r'))
source(get_path('code/plot.r'))
source(get_path('code/singlecell.r'))
source(get_path('code/scores.r'))
source(get_path('code/spatial.r'))
source(get_path('code/tpm.r'))
source(get_path('code/util.r'))
source(get_path('code/var_genes.r'))


# -------------
# load datasets
# -------------
# p = [spots x ident] cell frequencies
# x = [ident x genes] gene expression
# y = [spots x genes] gene expression (observed)
# e = [spots x genes] gene expression (predicted)

# visium data
sco = readRDS(get_path('sco.rds'))
vis = readRDS(get_path('vis.rds'))
gep = readRDS(get_path('1.BayesPrism/sco.gep.rds'))

# annotations
anno = ffread(get_path('anno_map.txt'), header=T, as.dt=T)
anno = anno[!is.na(finalrev)]
i2f = setNames(anno$finalrev, anno$annotation2rev)
f2g = setNames(anno$annotation3, anno$finalrev)

# fix cell type names
rownames(gep) = as.character(i2f[rownames(gep)])

# normalize gep matrix
gep.norm = 1e4*gep/rowSums(gep)

# deconv data
p = readRDS(get_path('1.BayesPrism/theta.rds'))
q = readRDS(get_path('1.BayesPrism/theta.merged.rds'))
x = readRDS(get_path('1.BayesPrism/all.Zkg.rds'))
X = readRDS(get_path('1.BayesPrism/all.Zkg.merged.rds'))
y = t(vis$counts[colnames(x),rownames(p)])
Y = t(vis$counts[colnames(X),rownames(q)])

# get samples
s = gsub('\\..*', '', rownames(y))

# scale data
x = 1e4*x/rowSums(x)
X = 1e4*X/rowSums(X)
y = 1e4*y/rowSums(y)

# predict expression
e = p %*% x
E = q %*% X
