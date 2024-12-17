# ----------------------------------------
# merge results across parallel lasso runs
# ----------------------------------------

b = c()
r = c()
v = c()

# load model data
# ---------------

for(i in 1:500){
    xi = readRDS(paste0('1.lasso_res.', i, '.500.rds'))
    xi = xi[lengths(xi) > 0]
    bi = sapply(xi, function(b) b$b)
    ri = sapply(xi, function(b) b$r)
    vi = sapply(xi, function(b) b$v.term/b$v.init)
    if(i == 1){
        b = bi
	r = ri
	v = vi
    } else {
        if(any(rownames(b) != rownames(bi))){stop()}
        if(any(rownames(r) != rownames(ri))){stop()}
        b = cbind(b, bi)
        r = cbind(r, ri)
        v = cbind(v, vi)
    }
}

saveRDS(b, file='beta.rds')
saveRDS(r, file='residuals.rds')
saveRDS(v, file='variance.rds')

