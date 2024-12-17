

mean_cv_loess = function(data, num_genes=1500, use_bins=TRUE, num_bins=20, window_size=100, do.plot=FALSE, invert=FALSE, use_var=FALSE){

    # calculate mean and cv
    u = apply(data, 1, mean)
    v = apply(data, 1, var)
    i = u > 0 & v > 0
    u = u[i]
    v = v[i]
    if(use_var == TRUE){
        cv = v
    } else {
        cv = sqrt(v)/u
    }

    # fit loess curve
    l = loess(log(cv) ~ log(u), family='symmetric')
    d = log(cv) - l$fitted

    # get variable genes
    if(use_bins == TRUE){

        # select variable genes from equal frequency bins
	k = as.integer(num_genes/num_bins)
	if(invert == FALSE){
	    var_genes = as.character(unlist(tapply(d, cut2(u, g=num_bins), function(a) names(sort(a, decreasing=T)[1:k]))))
	} else {
	    var_genes = as.character(unlist(tapply(d, cut2(u, g=num_bins), function(a) names(sort(a, decreasing=F)[1:k]))))
	}

    } else {

	# re-order by mean expression
	D = d[order(u)]

	# rolling z-scores (only use unique values)
	ru = rollapply(D, window_size, function(a) mean(unique(a)), partial=T)
	rs = rollapply(D, window_size, function(a) sd(unique(a)), partial=T)
	rz = structure((D - ru)/rs, names=names(D))

	# select variable genes
	var_genes = names(sort(rz, decreasing=T)[1:num_genes])
    }

    # plot results
    if(do.plot == TRUE){
	colors = c('#cccccc', 'black')[as.integer(names(u) %in% var_genes) + 1]
        plot(log(u), log(cv), pch=16, col=colors, xlab='log(mean)', ylab='log(cv)')
	lines(l$x[order(u)], l$fitted[order(u)], col='red', lw=2)
    }

    return(var_genes)
}


get_var_genes = function(data, ident=NULL, use_var=FALSE, genes.use=NULL, method='loess', num_genes=1500, min_cells=5, min_ident=25, do.plot=F, prefix=NULL, do.flatten=T, n.cores=1, ...){

    # Fix ident
    if(is.null(ident)){ident = rep(1, ncol(counts))}
    ident = as.factor(as.character(ident))

    # Filter ident
    u = sort(table(ident))
    i = names(which(u <= min_ident))
    if(sum(u[i]) >= min_ident){
        new_name = 'Merge'
    } else {
        new_name = names(which.min(u[u >= min_ident]))
    }
    levels(ident)[levels(ident) %in% i] = new_name
    print('Calculating var genes across')
    print(table(ident))
    
    # Subset idents
    levels.use = names(sort(table(ident), dec=T))[1:min(50, length(unique(ident)))]
    print(table(ident)[levels.use])

    # Subset genes
    if(is.null(genes.use)){
        genes.use = rownames(data)
    } else {
        genes.use = intersect(rownames(data), genes.use)
	print(paste('Subsetting from', nrow(data), 'to', length(genes.use), 'genes'))
    }
    data = data[genes.use,]

    # var_genes = sapply(levels(ident), function(i){print(i)
    var_genes = sapply(levels.use, function(i){print(i)

            # Start plotting device
	    if(!is.null(prefix)){png(paste(prefix, i, 'png', sep='.'), w=1000, h=800)}

	    # Subsample data
	    data = data[,ident == i]
    	    genes.use = rowSums(data > 0) >= min_cells
	    data = data[genes.use,]
	    print(dim(data))

	    if(method == 'loess'){
	        vi = mean_cv_loess(data, use_var=use_var, num_genes=num_genes, do.plot=do.plot, ...)
	    }
	    
	    # Stop plotting device
	    if(!is.null(prefix)){dev.off()}

	    return(vi)
    })

    if(do.flatten == TRUE){
        a = sort(table(as.character(unlist(var_genes))), decreasing=T)
	num_genes = min(length(a), num_genes)
	k = a[num_genes]
	u = names(a)[a > k]
	v = sample(names(a)[a == k], num_genes - length(u))
	var_genes = c(u,v)
    }

    return(var_genes)
}


