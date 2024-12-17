batch_correct = function(obj, batch, design=NULL, method='combat', genes.use='var', liger.data=NULL, ndim=10){

    # note: design matrix does not include batch
    if(!is.null(design)){
	if(is.null(rownames(design))){
            if(nrow(design) == ncol(obj$data)){
	        rownames(design) = colnames(obj$data)
	    } else {
	        stop('Error: check dimensions')
	    }
	}
        # subset design matrix
        design = design[colnames(obj$data),,drop=F]
    }

    # Get genes.use
    if(length(genes.use) == 1){
        if(genes.use == 'all'){
            print('Running batch correction on all genes')
	    genes.use = rownames(obj$data)
    	} else if(genes.use == 'var'){
            print('Running batch correction on variable genes')
	    if(length(obj$var.genes) == 0){
	        obj$var.genes = get_var_genes(obj, method='loess', num_genes=1500)
	    }
	    genes.use = obj$var.genes
    	} else {
            stop("Error: can only use 'all' or 'var' genes")
    	}
    }

    if(method == 'combat'){
	if(is.null(design)){
	    print('Running ComBat')
    	    print(table(batch))
	    new.data = ComBat(dat=as.matrix(obj$data[genes.use,]), batch, par.prior=T, prior.plots=F)
	} else {
	    print(sprintf('Fixing batches while controlling for %d covariates', ncol(design)))
	    print(table(batch))
	    model = model.matrix(~ ., data=design)
	    print(dim(model))
	    new.data = ComBat(dat=t(as.matrix(obj$data[genes.use,])), batch=batch, mod=model, par.prior=T, prior.plots=F)
	}
    }
        
    if(method == 'liger'){

	# Split data into batches
	sparse_split = function(x, g){
	    g = as.factor(g)
	    sapply(levels(g), function(gi){
	        x[g == gi,,drop=F]
	    }, simplify=F)
	}

	# Select batch
	batch = as.character(batch)
	batch[is.na(batch)] = 'Other'
	u = sort(table(batch))
	i = names(which(u <= 50))
	batch[batch %in% i] = 'Other'
	u = sort(table(batch))
	print(u)
	if(any(u <= 50)){
	    i = max(which(cumsum(u) <= 50)) + 1
	    i = names(u)[1:min(i, length(u))]
	    batch[batch %in% i] = 'Other'
	}
	batch = as.factor(batch)
	print(table(batch))

	new.data = obj$counts
	new.data = sapply(sparse_split(t(new.data), batch), t, simplify=F)
	for(i in 1:length(new.data)){
	    new.data[[i]][,1] = new.data[[i]][,1] + runif(nrow(new.data[[i]]), min=0, max=.01)
	}

	print('ndim')
	print(ndim)
	new.data = createLiger(new.data, remove.missing=F, take.gene.union=T)
	new.data = normalize(new.data)
	new.data@var.genes = obj$var.genes
	print(head(new.data@var.genes))
	new.data = scaleNotCenter(new.data)
	
	# Batch correct data
	new.data = optimizeALS(new.data, k=ndim)
	new.data = quantileAlignSNF(new.data)
	print(new.data@H.norm[1:5,1:5])
	print(head(colnames(obj$data)))
	print(table(colnames(obj$data) %in% rownames(new.data@H.norm)))
	new.data = new.data@H.norm[colnames(obj$data),]
    }

    return(as.data.frame(new.data))
}
