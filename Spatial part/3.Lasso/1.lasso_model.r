source('../load_data.r')
library(glmnet)

# ---------------
# input arguments
# ---------------
args = commandArgs(trailingOnly=T)
I = as.integer(args[[1]])
N = as.integer(args[[2]])
out = paste('1.lasso_res', I, N, 'rds', sep='.')
print(out)

# -------------
# load datasets
# -------------
# p = [spots x ident] cell frequencies
# x = [ident x genes] gene expression
# y = [spots x genes] gene expression (observed)
# e = [spots x genes] gene expression (predicted)

# -----------------------------
# get residuals from BayesPrism
# -----------------------------
e = readRDS(get_path('1.BayesPrism/total.Znkg.rds'))

# ----------------
# lasso regression
# ----------------

# pre-compute p %*% p
j = expand.grid(colnames(p), colnames(q))
P = p[,j[,1]] * q[,j[,2]]
colnames(P) = apply(j, 1, paste, collapse='.')
P = t(P)

# subset genes
chunk2 = function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))
genes.use = chunk2(colnames(y), N)[[I]]

res = sapply(genes.use, function(gi){tryCatch({
    
    # build features matrix [spots x ident^2]
    X = t(P * x[,gi])
    
    # initialize coefficients
    b = setNames(rep(0, ncol(X)), colnames(X))
    
    # select columns for lasso regression
    w = colSums(X)
    w = w/sum(w)
    j = names(which(w > .01))
    
    if(length(j) == 0){
        return(list())
    }
    
    # fit lasso regression
    fit = glmnet(X[,j,drop=F], y[,gi]-e[,gi], alpha=1, intercept=F, standardize=F, lower.limits=-1, upper.limits=100)
    b[j] = fit$beta[,ncol(fit$beta)]
    
    r = (y[,gi] - e[,gi]) - predict(fit, X[,j])[,ncol(fit$beta)]
    
    # calculate variances
    v.init = var(y[,gi] - e[,gi])
    v.final = var(r)
    v.term = (y[,gi] - e[,gi]) - t(t(X)*b)
    v.term = apply(v.term, 2, var)
    
    list(b=b, r=r, v.init=v.init, v.term=v.term, v.final=v.final)
    
}, error=function(e){list()})}, simplify=F)

saveRDS(res, file=out)
