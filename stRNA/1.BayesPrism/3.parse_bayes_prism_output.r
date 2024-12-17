source('../load_data.r')

# ---------------------------------------------------------------------------------------
# Parse the Bayes Prism output files to extract results for each cell type across samples
# Input = cell type ("all") or cell group ("merged")
# Output = theta, Zkg, Znkg files (see BayesPrism documentation)
# ---------------------------------------------------------------------------------------

# ---------------
# input arguments
# ---------------
args = commandArgs(trailingOnly=T)
ident.use = args[[1]]
type = args[[2]] # "all" or "merged"
if(! type %in% c('all', 'merged')){
    quit('error')
}
if(type == 'all'){
    suffix = 'Znkg'
}
if(type == 'merged'){
    suffix = 'Znkg.merged'
}
print(ident.use)
print(type)
print(suffix)

# ------------------------
# load the bayesprism data
# ------------------------
# u = all idents gene expression (sum over spots)
# x = ident gene expression [spots x genes]
# y = total gene expression [spots x genes]
# z = cell frequencies [spots x cell types]

fns = list.files(pattern='*ted.rds')
print(length(fns))

res = sapply(fns, function(fn){
    a = readRDS(fn)$res$first.gibbs.res
    if(type == 'all'){
        u = apply(a$Znkg, c(2,3), sum)
        x = a$Znkg[,ident.use,]
	y = apply(a$Znkg, c(1,3), sum)
    } else {
        u = apply(a$Znkg.merged, c(2,3), sum)
	x = a$Znkg.merged[,ident.use,]
	y = apply(a$Znkg.merged, c(1,3), sum)
    }
    z = a$gibbs.theta
    Z = a$theta.merged
    list(u=u, x=x, y=y, z=z, Z=Z)
}, simplify=F)

# ---------------------------------------
# aggregate gene signatures per cell type
# ---------------------------------------
u = as.matrix(sparse_matsum(sapply(res, function(a) a$u, simplify=F)))
x = as.matrix(sparse_rbind(sapply(res, function(a) a$x, simplify=F)))
y = as.matrix(sparse_rbind(sapply(res, function(a) a$y, simplify=F)))
z = as.matrix(sparse_rbind(sapply(res, function(a) a$z, simplify=F)))
Z = as.matrix(sparse_rbind(sapply(res, function(a) a$Z, simplify=F)))

# -------------------------------------
# save Bayes Prism results as RDS files
# -------------------------------------

saveRDS(u, file=paste('all', gsub('Znkg', 'Zkg', suffix), 'rds', sep='.'))
saveRDS(x, file=paste(ident.use, suffix, 'rds', sep='.'))
saveRDS(y, file=paste('total',   suffix, 'rds', sep='.'))

if(! file.exists('theta.rds')){saveRDS(z, file='theta.rds')}
if(! file.exists('theta.merged.rds')){saveRDS(Z, file='theta.merged.rds')}
