source('../load_data.r')
library(TED)

# ---------------
# read input data
# ---------------

# read input arguments
# --------------------
args = commandArgs(trailingOnly=T)
name = args[[1]] # cell type name
out = args[[2]] # output file

# load input data
# ---------------
# read the gep matrix
gep = readRDS(get_path('sco.gep.rds'))
# read the visium analysis object
vis = readRDS(get_path('vis.rds'))
print(dim(gep))

# load cell type names and annotations
# ------------------------------------
anno = read.table(get_path('anno_map.txt'), sep='\t', header=T)
fmap = ffread(get_path('final_names.txt'), row.names=1, header=T)
fmap = sapply(rownames(fmap), function(a) colnames(fmap)[which(fmap[a,] == 100)])

# --------------------------------------------------------
# filter cell types and genes prior to BayesPrism analysis
# --------------------------------------------------------

# subset cells
# ------------

# remove "contaminating" cells (e.g. doublets)
anno[grep('^Contam', anno$final),][] = NA 

# select cells from given cell type
x = t(vis$counts[,grep(name, colnames(vis$counts))])
print(dim(x))

# subset genes
# ------------

# select the 2500 most highly expressed genes
g1 = names(sort(apply(gep, 2, max), dec=T)[1:2500])

# select the 5000 most highly variable genes
g2 = get_var_genes(t(gep), rep('a', nrow(gep)), num_genes=5000)

# select highly expressed and highly variable genes
genes.use = unique(sort(c(g1,g2)))
genes.use = intersect(intersect(genes.use, colnames(gep)), colnames(x))
gep = gep[,genes.use]
x = x[,genes.use]

# --------------
# run BayesPrism
# --------------

# get cell type labels
ident.use = setNames(anno$annotation3, anno$annotation2rev)[rownames(gep)]

# get cell subtype labels
rownames(gep) = as.character(fmap[rownames(gep)])

# run BayesPrism
res = run.Ted(ref.dat=gep,
              X = x,
	      cell.type.labels=ident.use,
	      cell.subtype.labels=rownames(gep),
	      input.type='GEP',
	      n.cores=16,
	      pdf.name=name)

# save BayesPrism output as RDS file
saveRDS(res, file=out)
