options(max.print=1000)

msg = function(text, verbose){
    if(verbose == TRUE){
        print(text)
    }
}

set.ident = function(obj, ident.use){
    obj$ident = setNames(as.factor(ident.use), colnames(obj$data))
    return(obj)
}

init_obj = function(){
    obj = list()
    obj$counts = NA
    obj$data = NA
    obj$meta.data = NA
    obj$ident = NA
    obj$tsne.rot = NA
    obj$image = NA
    obj$pca.obj = NA
    obj$pca.rot = NA
    return(obj)
}

make_obj = function(counts=NULL, regex='', regexv='', minc=10, maxc=NULL,maxc_per_group=NULL, ming=500, maxg=1e6, genes.use=NULL, cells.use=NULL, ident_fxn=NULL, verbose=FALSE, x11=FALSE, bayes.n=0, qnorm=F){

    # Load packages
    
    # Set graphics device
    options(device=pdf)

    # Load counts from singlecell object, matrix, or file
    msg('loading counts', verbose)
    if(typeof(counts) == typeof('')){
        if(file.exists(counts)){
	    counts = fread(paste('zcat', counts))
	    counts = data.frame(counts, row.names=1)
	} else {
	    counts = read_mtx(prefix=counts)
	}
    }
    msg( sprintf('Counts = %d x %d', nrow(counts), ncol(counts)), verbose)

    # Subset counts with NA, regex, genes.use, and cells.use
    msg( 'Subsetting counts', verbose)
    if(regex != ''){
        j = grep(regex, colnames(counts))
	counts = counts[,j,drop=F]
    }
    if(regexv != ''){
        j = grep(regexv, colnames(counts), invert=T)
	counts = counts[,j,drop=F]
    }
    if(!is.null(genes.use)){
	genes.use = intersect(rownames(counts), genes.use)
	counts = counts[genes.use,,drop=F]
    }
    if(!is.null(cells.use)){
	cells.use = intersect(colnames(counts), cells.use)
	counts = counts[,cells.use,drop=F]
    }
    genes.use = rowSums(is.na(counts)) == 0
    counts = counts[genes.use,,drop=F]

    msg( sprintf('counts = %d x %d', nrow(counts), ncol(counts)), verbose)
    
    # Convert counts to sparse matrix
    if(is.data.frame(counts)){counts = as.matrix(counts)}
    counts = as(counts, 'sparseMatrix')
    
    # Get cell identities
    if(is.null(ident_fxn)){
        ident = sapply(strsplit(colnames(counts), '\\.'), '[[', 1)
    } else {
        ident = sapply(colnames(counts), ident_fxn)
    }
    ident = setNames(as.factor(ident), colnames(counts))

    # Downsample cells
    if(!is.null(maxc) | !is.null(maxc_per_group)){
        cells.use = simple_downsample(cells=colnames(counts), groups=ident, total_cells=maxc, cells_per_group=maxc_per_group)
	msg( paste('Downsampling to', length(cells.use), 'total cells'), verbose)
	counts = counts[,cells.use]
    }

    # Filter cells by minc, ming, maxg
    msg( 'Filtering counts', verbose)
    j1 = colSums(counts > 0) >= ming
    j2 = colSums(counts > 0) <= maxg
    counts = counts[,(j1 & j2)]
    i = rowSums(counts > 0) >= minc
    counts = counts[i,]
    msg( sprintf('counts = %d x %d', nrow(counts), ncol(counts)), verbose)

    # Add Bayesian prior
    if(bayes.n > 0){
        print('Bayesian prior')
        p = scaleMargins(counts, cols=1/colSums(counts))
	p = rowMeans(p)
	p = rmultinom(n=ncol(counts), size=bayes.n, prob=p)
	print(paste('Adding', mean(colSums(p)), 'reads to every cell'))
	counts = counts + p
    }
    
    # Make singlecell object
    msg( 'Making singlecell object', verbose)
    obj = init_obj()
    obj$counts = counts
    obj$ident = ident[colnames(counts)]
    obj$meta.data = data.frame(nGene=colSums(counts > 0), nUMI=colSums(counts))
    rownames(obj$meta.data) = colnames(counts)
    
    # Normalize data (default = TPM)
    msg( 'Normalizing data', verbose)
    if(qnorm == FALSE){
        obj$data = calc_tpm(counts=obj$counts)
	obj$data@x = log2(obj$data@x + 1)
    }
        
    # Get cell identities
    msg( 'Setting cell identities', verbose)
    if(!is.null(ident_fxn)){
        ident = sapply(colnames(obj$data), ident_fxn)
	obj = set.ident(obj, ident.use=ident)
	obj$meta.data$orig.ident = obj$ident
    }

    if(length(unique(obj$ident)) > 100){
        msg( 'WARNING: nlevels(obj$ident) > 100', verbose)
        #obj$ident = '1'
    	#obj$meta.data$orig.ident = '1'
    }

    print(table(obj$ident))
    return(obj)
}

run_seurat = function(name, obj=NULL, counts=NULL, regex='', regexv='', cells.use=NULL, genes.use=NULL, minc=5, maxc=1e6, ming=200, maxg=1e6, ident_fxn=NULL, varmet='loess', var_regexv=NULL,
             var_remove=NULL, min_cv2=.25, var_genes=NULL, use_var=FALSE, var_add=NULL, bayes.n=0, qnorm=F, num_genes=1500, regress=NULL, do.batch='none', batch.use=NULL, store.batch=FALSE, bc.data=NULL,
	     design=NULL, pc.data=NULL, num_pcs=0, pcs.use=NULL, pcs.rmv=NULL, robust_pca=F, scale.max=10, 
	     perplexity=25, max_iter=1000, dist.use='euclidean', do.largevis=FALSE, do.umap=FALSE, largevis.k=50, do.fitsne=FALSE, fitsne.K=-1,
	     cluster='infomap', k=c(), verbose=T, write_out=T, do.backup=F, ncores=1, stop_cells=50, marker.test=''){

    # check input arguments
    if(! do.batch %in% c('none', 'combat', 'mnn', 'cca', 'multicca', 'liger', 'harmony')){stop('do.batch must be none, combat, mnn, or cca')}
    if(! is.null(batch.use)){if(is.null(names(batch.use))){stop('batch.use needs names')}}

    # Make singlecell object
    if(is.null(obj)){
        obj = make_obj(counts=counts, regex=regex, regexv=regexv, minc=minc, maxc=maxc, ming=ming, maxg=maxg, genes.use=genes.use, cells.use=cells.use, ident_fxn=ident_fxn, verbose=verbose, qnorm=qnorm, bayes.n=bayes.n)
    }
    if(ncol(obj$data) <= stop_cells){return(obj)}
    
    msg( 'Selecting variable genes', verbose)
    ident = obj$ident
    if(is.null(var_genes)){
        gi = rownames(obj$counts)
	print(paste('Starting with', length(gi), 'genes'))
	if(!is.null(var_regexv)){gi = grep(var_regexv, gi, invert=T, value=T)}
	print(paste('var_regexv:', length(gi), 'genes'))
	if(!is.null(var_remove)){gi = setdiff(gi, var_remove)}
	print(paste('var_remove:', length(gi), 'genes'))
	var_genes = get_var_genes(obj$counts, ident=ident, method=varmet, genes.use=gi, num_genes=num_genes, min_ident=25, use_var=use_var)
	if(!is.null(var_add)){
	    var_genes = c(var_genes, setdiff(intersect(gi, var_add), var_genes))
	}
    }
    
    #if(is.null(var_genes)){var_genes = get_var_genes(obj$counts, ident=ident, method=varmet, num_genes=num_genes, min_ident=25)}
    #if(!is.null(var_regexv)){var_genes = grep(var_regexv, var_genes, invert=T, value=T)}
    msg( sprintf('Found %d variable genes', length(var_genes)), verbose)
    obj$var.genes = intersect(var_genes, rownames(obj$data))
    print(var_genes)
    
    # Regression
    if(!is.null(regress)){
        print('Regressing out gene signature from pc.data')
	print(obj$data[1:5,1:5])
        x = score_cells(obj, names=regress)
	obj$data = t(apply(obj$data, 1, function(a) .lm.fit(as.matrix(x), as.matrix(a))$residuals[,1]))
	print(obj$data[1:5,1:5])
    }

    if(is.null(batch.use)){
	batch.use = obj$ident
    }
    batch.use = batch.use[names(obj$ident)]
    print(table(batch.use))
    
    # Batch correction with variable genes
    if(! do.batch %in% c('none', 'harmony')){
	msg( 'Batch correction', verbose)
	if(!is.null(design)){
	    if(length(intersect(names(obj$ident), rownames(design))) < 10){
	        rownames(design) = names(obj$ident)
	    }
	    design = design[names(obj$ident),,drop=F]
	}
	if(is.null(bc.data)){
	    bc.data = batch_correct(obj, batch.use, design=design, method=do.batch, genes.use=obj$var.genes, ndim=num_pcs)
	}
	print(dim(bc.data))	
	# store batch corrected data
	if(store.batch == TRUE){
	    obj$bc.data = bc.data
	}
	
	# write batch corrected data to file
	if(write_out == TRUE){fwrite(as.data.table(bc.data), file=paste0(name, '.bc.data.txt'), sep='\t')}
	
	pc.data = t(scale(t(bc.data), center=F))
	#pc.data = bc.data
    }
    if(is.null(pc.data)){pc.data = obj$data}

    # Number of significant PCs (stored in obj$meta.data$num_pcs)
    if(num_pcs == 0){
	num_pcs = sig.pcs.perm(scale(t(obj$data)), randomized=T, n.cores=ncores)$r + 2
    }
    if(is.na(num_pcs)){num_pcs = 5}
    obj$meta.data$num_pcs = num_pcs
    msg( sprintf('Found %d significant PCs', num_pcs), verbose)
    
    # Fast PCA on data
    if(do.batch %in% c('liger', 'multicca')){
        obj$pca.rot = as.data.frame(pc.data) # liger and multi-cca output saved in pc.data
	print(dim(obj$pca.rot))
    } else {
        obj = run_rpca(obj, data=pc.data, k=50, genes.use=obj$var.genes, robust=robust_pca, rescale=T, scale.max=scale.max)
	if(write_out == TRUE){saveRDS(obj$pca.obj, file=paste0(name, '.pca.rds'))}
    }
    
    # Harmony batch correction
    if(do.batch == 'harmony'){
	print('Running harmony')
	print(table(batch.use))
	obj$pca.rot = HarmonyMatrix(obj$pca.rot[,1:num_pcs], batch.use, 'dataset', do_pca=FALSE)
    }
        
    # Fix problem with duplicates
    obj$pca.rot[,num_pcs] = obj$pca.rot[,num_pcs] + runif(nrow(obj$pca.rot), min=-1e-8, max=1e-8)

    # Regress out PCs
    if(is.null(pcs.use)){
        pcs.use = 1:(min(ncol(obj$pca.rot), num_pcs))
    }

    # TSNE
    knn = NULL
    if(max_iter > 0){
    if(do.fitsne == TRUE){
        msg( 'FIt-SNE', verbose)
	q = fftRtsne(obj$pca.rot[,pcs.use], max_iter=max_iter, perplexity=perplexity, K=fitsne.K, fast_tsne_path='~/code/FIt-SNE/bin/fast_tsne')
    } else if(do.largevis == TRUE){
        msg( 'largeVis', verbose)
	q = run_largevis(obj$pca.rot[,pcs.use], k=largevis.k, save_knn=TRUE, save_weights=FALSE, dist.use=dist.use, verbose=T)
	knn = q$knn
	q = q$coords
    } else if(do.umap == TRUE){
        msg( 'umap', verbose)
	q = umap(obj$pca.rot[,pcs.use], method='umap-learn')$layout
	print(dim(q))
	#q = umap(t(as.matrix(pc.data)), method='umap-learn')$layout
    } else {
        msg( 'TSNE', verbose)
        if(dist.use == 'euclidean'){
            q = Rtsne(obj$pca.rot[,pcs.use], do.fast=T, max_iter=max_iter, perplexity=perplexity, verbose=T)$Y
        } else if(dist.use == 'cosine'){
            d = cosine_dist(t(obj$pca.rot[,pcs.use]))
	    q = Rtsne(d, is_distance=T, do.fast=T, max_iter=max_iter, perplexity=perplexity, verbose=T)$Y
        } else {
            d = dist(obj$pca.rot[,pcs.use], method=dist.use)
	    q = Rtsne(d, is_distance=T, do.fast=T, max_iter=max_iter, perplexity=perplexity, verbose=T)$Y
        }
    }
    rownames(q) = colnames(obj$data)
    colnames(q) = c('tSNE_1', 'tSNE_2')
    obj$tsne.rot = as.data.frame(q)
    }

    # Cluster cells and run DE tests
    if(length(k) > 0){

        # Save backup singlecell object
	if(do.backup){saveRDS(obj, file=paste0(name, '.obj.rds'))}

	msg( 'Clustering cells', verbose)
	k = k[k < ncol(obj$data)]
	u = paste('Cluster.Infomap.', k, sep='')
	v = run_graph_cluster(data=obj$pca.rot[,pcs.use], k=k)
	obj$meta.data[,u] = v

	if(marker.test != ''){
    	    msg( 'Differential expression', verbose)
            covariates = subset(obj$meta.data, select=c(nGene, Cell_Cycle))
	    obj = set.ident(obj, ident.use=v[,1])
	    print(table(obj$ident))
	    markers = lapply(k, function(ki){
	        obj = set.ident(obj, ident.use=obj$meta.data[,paste0('Cluster.Infomap.', ki)])
	        markers = p.find_all_markers(obj, test.use=marker.test)
	    })
	    names(markers) = k
	}
    }

    if(write_out){

	# Plot TSNE
	png(paste0(name, '.tsne.png'), width=800, height=650)
	plot_tsne(obj, pt.size=1)
	dev.off()

	# Plot clusters
	if(length(k) > 0){
	    pdf(paste0(name, '.clusters.pdf'), width=9, height=9)
	    plot_clusters(obj)
	    dev.off()
	}

	# Marker genes
	if(marker.test != ''){
	    for(ki in names(markers)){write.table(markers[[ki]], file=paste0(name, '.k', ki, '.', marker.test, '.txt'), sep='\t', quote=F)}
	}

	# Save singlecell object
	saveRDS(obj, file=paste0(name, '.obj.rds'))
    }

    return(obj)
}
