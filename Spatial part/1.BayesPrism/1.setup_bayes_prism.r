source('../load_data.r')

# -----------------------------------------------------------------
# Calculate a gene expression program (gep) matrix, which describes
# for each cell type (rows) the mean expression of all genes (cols)
#
# The output is an RDS file which contains the gep matrix:
# gep = [cell types x genes]
# -----------------------------------------------------------------

# Load "excluded" genes from BayesPrism directory
glist = read.table('genelist.hs.txt', sep='\t')
glist.rm = sort(unique(as.character(glist[,3])))

# ---------------------------------------------------------
# Calculate the mean expression profile for every cell type
# ---------------------------------------------------------

# Load the single-cell analysis object
sco = readRDS(get_path('sco.rds'))

# Filter gene set:
# ----------------
# Only use genes found in > 100 cells
genes.use = rownames(sco$counts)[rowSums(sco$counts > 0) >= 100]

# Exclude genes in the BayesPrism gene set
genes.use = setdiff(intersect(rownames(sco$counts), genes.use), glist.rm)

# Exclude "contaminating" cell types (doublets, etc)
cells.use = names(sco$ident)[!grepl('Contam|Ribo', sco$ident)]

# Calculate the mean expression profiles for each cell type
gep = nice_agg2(t(sco$counts[genes.use,cells.use]), droplevels(sco$ident[cells.use]), mean)

# Save GEP matrix as RDS file
saveRDS(gep, file='sco.gep.rds')
