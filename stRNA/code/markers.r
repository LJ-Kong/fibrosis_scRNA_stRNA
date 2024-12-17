
load_mast = function(){
    library(BiocGenerics, pos=length(search()))
    library(S4Vectors, pos=length(search()))
    library(DelayedArray, pos=length(search()))
    library(MAST)
}


test_log_base = function(obj, base=2, total=1e4){
    # Test log base of obj$data
    j = sample(1:ncol(obj$data), 1)
    u = sum(base**obj$data[,j] - 1)
    (1e4 - 1e-2) <= u & u <= (1e4 + 1e-2)
}


unorder_factors = function(x){
    j = sapply(x, is.ordered)
    x[,j] = lapply(x[,j,drop=F], function(a) factor(as.character(a), levels=levels(a)))
    x
}

relevel_factors = function(x){
    j = sapply(x, is.factor)
    x[,j] = lapply(x[,j,drop=F], function(a) droplevels(a))
    j = sapply(x, function(a) length(unique(a)) > 1)
    x = x[,j,drop=F]
    x
}


get_data = function(obj, data.use='tpm', tpm=NULL, cells.use=NULL){

    # Retrieve data from a singlecell object
    # data.use can be: counts, tpm, log2, data, any matrix
    # can pre-calculate tpm for speed
    
    if(!is.character(data.use)){
        return(data.use)
    }

    if(data.use == 'counts'){
        data = obj$counts
    }

    if(data.use == 'tpm'){
        if(is.null(tpm)){
	    print(dim(data))
	    data = calc_tpm(obj, cells.use=cells.use)
	    print(dim(data))
	} else {
	    data = tpm
	}
    }
    
    if(data.use == 'log2'){
        if(!test_log_base(obj, base=2, total=1e4)){
	    stop('Error: obj$data log base != 2')
        }
        data = obj$data
    }

    return(data)
}


select_cells = function(obj, covariates, batch.use=NULL, cells.use=NULL, max_cells=NULL){

    # Select cells from groups defined by covariates matrix
    # -----------------------------------------------------
    # 1. select cells.use
    # 2. build groups from covariates matrix
    # 3. select max_cells from each group
    # 4. sample evenly across batch.use

    cat('\nSelecting cells\n')
    
    # Get batch
    ncells = ifelse(is.null(ncol(obj$data)), nrow(covariates), ncol(obj$data))
    if(is.null(ncells)){ncells = nrow(covariates)}
    if(is.null(batch.use)){batch.use = rep('All', ncells)}
    
    # Subset data
    i = apply(covariates, 1, function(a) !any(is.na(a)))
    covariates = covariates[i,,drop=F]
    batch.use = batch.use[i]
    
    # Get cells to use
    cells = rownames(covariates)
    if(!is.null(cells.use)){
        cells = intersect(cells, cells.use)
    }
    
    # Construct cell groups
    j = sapply(covariates, function(a) !is.numeric(a))
    groups = as.factor(apply(covariates[cells,j,drop=F], 1, function(a) paste(a, collapse='.')))
    cat('\n')
    print(table(groups))

    # Select max_cells from each group
    if(!is.null(max_cells)){

        # Combine batch.use and groups to sample evenly across batches
        batch.use = as.factor(paste0(groups, batch.use))
	print(table(batch.use))

	# Sample [max_cells] separately from each group
	cells = unname(unlist(lapply(levels(groups), function(group){
	    simple_downsample(cells=cells[groups == group], groups=batch.use[groups == group], total_cells=max_cells)
	})))

	# Print subsampled group sizes
	groups = apply(covariates[cells,j,drop=F], 1, function(a) paste(a, collapse='.'))
	cat('\n')
	print(table(groups))
    }
    return(cells)
}


select_genes = function(obj, stats, data.use=NULL, genes.use=NULL, min_cells=3, min_alpha=.025, min_fc=1.2, dir='both'){

    # Select genes from singlecell object by cells, genes.use, min_cells, min_alpha, and min_fc

    # Select genes by genes.use
    if(is.null(data.use)){data.use = obj$data}
    if(is.null(genes.use)){genes.use = rownames(data.use)}
    
    # Check stats
    if(is.null(stats)){
        
	g1 = names(which(rowSums(obj$counts > 0) >= min_cells))
	g2 = names(which(rowMeans(obj$counts > 0) >= min_alpha))
	genes.use = intersect(g1, g2)
	
    } else {
        
    	# Select genes by min_cells
    	g1 = as.character(stats[, (max(n) >= min_cells) | (max(ref_n) >= min_cells), .(gene)][V1 == TRUE, gene])
	
    	# Select genes by min_alpha
    	g2 = as.character(stats[, (max(alpha) >= min_alpha) | (max(ref_alpha) >= min_alpha), .(gene)][V1 == TRUE, gene])
	
    	# Select genes by min_fc
    	if(dir == 'pos'){
            g3 = as.character(stats[, max(log2fc) >= log2(min_fc), .(gene)][V1 == TRUE, gene])
    	} else {
            g3 = as.character(stats[, max(abs(log2fc)) >= log2(min_fc), .(gene)][V1 == TRUE, gene])
        }
	
        # Intersect and return genes.use
    	genes.use = Reduce(intersect, list(genes.use, g1, g2, g3))
    }
    
    return(genes.use)
}


expression_stats = function(tpm, covariates, formula, lrt_regex, genes.use=NULL, cells.use=NULL, invert_method='auto', invert_logic='last'){

    # Calculate expression statistics for groups specified by formula
    # each column of the model matrix is either a single term A or an interaction term A:B
    # if single term, then select cells A and ~A
    # if interaction, then select cells A:B and A:~B
    # for multi-level factors, ~A is the next highest level of A
        
    # Select cells
    if(!is.null(cells.use)){
        tpm = tpm[,cells.use]
	covariates = covariates[cells.use,,drop=F]
    }
    print(sum(tpm))
    
    # Drop levels
    covariates = relevel_factors(covariates)
    
    # Check for factors
    u = grepl(lrt_regex, colnames(covariates))
    v = sapply(covariates, is.factor)
    
    if(sum(u & v) == 0){
        print('Skipping expression statistics')
        return(NULL)
    }
       
    # Select genes
    if(!is.null(genes.use)){tpm = tpm[genes.use,,drop=F]}
    total = rowSums(tpm)
    
    # Model matrix
    print('Model matrix')
    mm_formula = gsub('\\+ *\\S*\\|\\S+', '', as.character(formula), perl=T)
    mm = as.matrix(model.matrix(as.formula(mm_formula), data=unorder_factors(covariates)))
    
    # Invert matrix
    print('Invert matrix')
    u = mm_logical_not(mm, formula, covariates, method=invert_method, invert=invert_logic)
    MM = u$x
    refs = structure(u$names, names=colnames(MM))

    # For every column that matches lrt_regex
    print('Expression stats')
    stats = lapply(grep(lrt_regex, colnames(mm), value=T), function(a){print(a)
        
	# cell indices
	i = as.logical(mm[,a])
	j = as.logical(MM[,a])
	ref = refs[[a]]

	# number of expressing cells
	n1 = rowSums(tpm[,i,drop=F] > 0)
	n2 = rowSums(tpm[,j,drop=F] > 0)

	# total expression over all cells
	s1 = rowSums(tpm[,i,drop=F])
	s2 = rowSums(tpm[,j,drop=F])

	# fraction of expressing cells (alpha)
	a1 = n1/sum(i)
	a2 = n2/sum(j)

	# mean over expressing cells (mu)
	m1 = ifelse(n1 > 0, s1/n1, 0)
	m2 = ifelse(n2 > 0, s2/n2, 0)

	# mean over all cells
	u1 = s1/sum(i)
	u2 = s2/sum(j)

	# fraction of total expression
	t1 = s1/total
	t2 = s2/total

	# fix zeros for logs
	if(class(tpm) == 'dgCMatrix'){
	    zero = .5*min(tpm@x)/(sum(i) + sum(j)) # fast min (sparse matrix)
	} else {
	    zero = .5*min(tpm)/(sum(i) + sum(j)) # slow min (other matrix)
	}
	m1 = m1 + .5*zero
	m2 = m2 + .5*zero
	u1 = u1 + .5*zero
	u2 = u2 + .5*zero

	# log fold change
	log2fc = log2(u1) - log2(u2)

	# combine in data frame
	res = data.frame(gene=rownames(tpm), contrast=a, ref=ref, n=n1, ref_n=n2, alpha=a1, ref_alpha=a2, mu=log2(m1), ref_mu=log2(m2), mean=log2(u1), ref_mean=log2(u2), total=t1, ref_total=t2, log2fc=log2fc)
	return(res)
    })
    stats = as.data.table(do.call(rbind, stats))
    return(stats)
}


fdr_stats = function(data, covariates, formula, lrt_regex, genes.use=NULL, cells.use=NULL, invert_method='multi', invert_logic='last'){

    # Calculate FDR statistics for groups specified by formula
    
    # Select genes
    if(is.null(genes.use)){genes.use = rownames(data)}
    if(is.null(cells.use)){cells.use = colnames(data)}
    data = data[genes.use, cells.use]
    covariates = covariates[cells.use,,drop=F]

    # Check for factors
    u = grepl(lrt_regex, colnames(covariates))
    v = apply(covariates, 2, is.factor)
    if(sum(u & v) == 0){
        print('Skipping FDR statistics')
        return(NULL)
    }
        
    # Model matrix
    mm = as.matrix(model.matrix(as.formula(formula), data=unorder_factors(covariates)))

    # Invert matrix
    u = mm_logical_not(mm, formula, covariates, method=invert_method, invert=invert_logic)
    MM = u$x
    refs = structure(u$names, names=colnames(MM))

    # For every column that matches lrt_regex
    stats = lapply(grep(lrt_regex, colnames(mm), value=T), function(a){

	# Select cells in each group
	i = as.logical(mm[,a])
	j = as.logical(MM[,a])
	ref = refs[[a]]

	# False discovery rate
	fdr = t(apply(data, 1, function(b){
	    predictions = b[i|j]
	    labels = ifelse(i, TRUE, NA)
	    labels[j] = FALSE
	    labels = na.omit(labels)
	    calc_fdr(predictions, factor(labels, levels=c(FALSE, TRUE)))
 	}))
	colnames(fdr) = c('cutoff', 'accuracy', 'sensitivity', 'specificity', 'fdr', 'f1')
	fdr = data.frame(gene=rownames(data), contrast=a, ref=ref, fdr)
	return(fdr)
    })
    stats = as.data.table(do.call(rbind, stats))
    return(stats)
}


p.find_markers = function(obj=NULL, ident.1=NULL, ident.2=NULL, ident.use=NULL, tpm.use='tpm', data.use='log2', genes.use=NULL, cells.use=NULL, test.use='roc', min_cells=3, min_alpha=.05, min_fc=1.25,
                          max_cells=1000, batch.use=NULL,
	       	          dir='pos', tpm=NULL, covariates=NULL, formula='~ ident', lrt_regex='ident', gsea.boot=100, invert_method='auto', invert_logic='last', do.stats=FALSE, n.cores=1,
			  filter_genes=TRUE, qnorm=FALSE, glmer=FALSE, fix_rank=FALSE, skip_tpm=FALSE){

    # Allow running without obj
    if(is.null(obj)){
        if(is.null(nrow(data.use))){quit('must specify obj or data.use')}
        tpm.use = data.use
    }

    # Get cell identities
    if(is.null(ident.use)){ident.use = obj$ident}
    
    # Fix cells.use to allow boolean
    if(typeof(cells.use) == 'logical'){
        cells.use = colnames(obj$data)[cells.use]
    }
    
    # Check lrt regex
    if(!grepl(lrt_regex, formula)){
        print('lrt_regex does not match formula')
	return(FALSE)
    }
    
    print(head(ident.use))
    # Build covariates
    if(!is.null(ident.1)){
	if(is.null(ident.2)){
	    ident.use = as.factor(ifelse(ident.use == ident.1, ident.1, 'Other'))
	    ident.use = relevel(ident.use, 'Other')
	} else {
	    ident.use = ifelse(ident.use %in% c(ident.1, ident.2), as.character(ident.use), NA)
	    ident.use = factor(ifelse(ident.use %in% c(ident.1, ident.2), as.character(ident.use), NA), levels=c(ident.2, ident.1))
	}
	if(is.null(covariates)){
	    covariates = data.frame(ident=ident.use)
	} else {
	    covariates$ident = ident.use
	}
    }
    if(!is.null(obj)){
        rownames(covariates) = colnames(obj$data)
    } else {
        rownames(covariates) = colnames(data.use)
    }
    print(table(covariates$ident))
    print(dim(covariates))
    print(head(covariates))
    
    # Check covariates
    q = sapply(covariates, typeof)
    if('character' %in% q){print(q); stop('error: invalid covariates type')}
    
    # Select cells
    print('selecting cells')
    cells.use = select_cells(obj, covariates, cells.use=cells.use, max_cells=max_cells, batch.use=batch.use)
    
    # TPM for log fold changes [genes x cells]
    print('getting data')
    if(!skip_tpm){
        tpm = get_data(obj, data.use=tpm.use, tpm=tpm, cells.use=cells.use)
    }
    
    # Data for DE test [genes x cells]
    data = get_data(obj, data.use=data.use, tpm=tpm, cells.use=cells.use)
    
    print('Expression stats')
    if(skip_tpm){
        print('Skipping TPM... expression_stats will be incorrect')
        tpm = data
    }    
    stats = expression_stats(tpm, covariates, formula, lrt_regex, genes.use=genes.use, cells.use=cells.use, invert_method=invert_method, invert_logic=invert_logic)
    print(dim(stats))
    
    if(filter_genes == TRUE){

        # Expression stats
        #stats = expression_stats(tpm, covariates, formula, lrt_regex, genes.use=genes.use, cells.use=cells.use, invert_method=invert_method, invert_logic=invert_logic)

	# Select genes
        genes.use = select_genes(obj, stats, data.use=data, genes.use=genes.use, min_cells=min_cells, min_alpha=min_alpha, min_fc=min_fc, dir=dir)
        if(length(genes.use) == 0){return(c())}

    } else {#stats = NULL;
        if(is.null(genes.use)){genes.use = rownames(data)}
    }
    print(length(genes.use))

    # FDR stats
    if(do.stats == TRUE){
        print('FDR stats')
        fdr = fdr_stats(data, covariates, formula, lrt_regex, genes.use=genes.use, cells.use=cells.use, invert_method=invert_method, invert_logic=invert_logic)
	setkey(stats, gene, contrast, ref)
	setkey(fdr, gene, contrast, ref)
	stats = stats[fdr,]
    }

    # Subset data
    print(paste('Testing', length(genes.use), 'genes in', length(cells.use), 'cells'))
    data = data[genes.use, cells.use, drop=F]
    data = rbind(data, rnorm(ncol(data)))
    covariates = unorder_factors(covariates[cells.use, , drop=F])
    covariates = relevel_factors(covariates)

    # Fix rank deficiency (zero lowest expressing cell for each gene)
    if(fix_rank == TRUE){
        for(gi in rownames(data)){
	    tryCatch({
	        j = data[gi,] > 0
	        u = names(which.min(data[gi, j]))
	        v = names(which.max(data[gi, !j]))
	        data[gi, v] = data[gi, u]
	        data[gi, u] = 0
	    }, error=function(e){NULL})
	}
    }

    print(dim(data))

    # Run marker tests
    labels = covariates[,1]
    if(test.use == 'f'){markers = de.rocr(data, labels, measures='f')}
    if(test.use == 'fdr'){markers = de.fdr(data, labels)}
    if(test.use == 'pr'){markers = de.rocr(data, labels, measures=c('prec', 'rec'))}
    if(test.use == 'mast'){markers = de.mast(data, covariates=covariates, formula=formula, lrt_regex=lrt_regex, n.cores=n.cores, glmer=FALSE)}
    if(test.use == 'gsea'){markers = gsea.mast(data, covariates=covariates, formula=formula, lrt_regex=lrt_regex, gsea.boot=gsea.boot, n.cores=n.cores)}
    if(test.use == 'roc'){markers = de.rocr(data, labels, measures='auc')}
    if(test.use == 'mpath'){markers = de.mpath(obj$counts[genes.use,cells.use], covariates=covariates, formula=formula)}

    # Add cluster information
    if(! 'contrast' %in% colnames(markers)){
        markers$contrast = paste0('ident', levels(labels)[nlevels(labels)])
    }

    # Merge results
    markers = as.data.table(markers)
    setkey(markers, gene, contrast)
    if(!is.null(stats)){
        setkey(stats, gene, contrast)
        markers = markers[stats,]
    }

    # Sort markers
    if('auc' %in% colnames(markers)){
        markers = markers[order(contrast, -1*auc),]
    } else if('f1' %in% colnames(markers)){
        markers = markers[order(contrast, -1*f1),]
    } else if('pval' %in% colnames(markers)){
        markers = markers[order(contrast, pval),]
    } else if('pvalH' %in% colnames(markers)){
        markers = markers[order(contrast, pvalH),]
    } else {
        markers = markers[order(contrast, -1*alpha),]
    }

    # Add ident
    markers$ident = paste(ident.1, ident.2, sep=';')
    markers$ident = gsub(';$', '', markers$ident)

    # Return marker genes
    return(markers)
}

p.find_all_markers = function(obj, idents=NULL, ident.use=NULL, data.use='log2', tpm=NULL, do.precalc=T, n.cores=1, qnorm=FALSE, ...){

    # note: idents = list of idents to test (default = all)
    
    # Get cell identities
    if(is.null(ident.use)){ident.use = obj$ident}
    ident.use = as.factor(ident.use)
    if(is.null(idents)){idents = as.character(levels(ident.use))} else {idents = unique(as.character(idents))}
    print(paste('Calculating markers for:', paste(idents, collapse=', ')))

    # Pre-calculate TPM and data
    if(do.precalc == TRUE){
        tpm = get_data(obj, data.use='tpm', tpm=tpm)
        data.use = get_data(obj, data.use=data.use, tpm=tpm)
    }
    
    # Find marker genes
    do.call(rbind, sapply(idents, function(i) p.find_markers(obj, ident.1=i, ident.use=ident.use, tpm=tpm, data.use=data.use, n.cores=1, qnorm=FALSE, ...), simplify=F))

}

calc_rocr = function(predictions, labels, measures='auc', retx=FALSE){
    
    if(length(measures) == 1){
        q = performance(prediction(predictions, labels, label.ordering=levels(labels)), measures)
    } else {
        q = performance(prediction(predictions, labels, label.ordering=levels(labels)), measures[[1]], measures[[2]])
    }
    
    x = q@x.values
    if(length(x) > 0){x = ifelse(is.na(x[[1]]), 0, x[[1]])}
    y = ifelse(is.na(q@y.values[[1]]), 0, q@y.values[[1]])
    
    if(length(y) > 1){
        if(q@x.name == 'Cutoff'){
	    i = which.max(y[x != min(x)])
	    x = x[i]
	    y = y[i]
	} else {
            y = sum(diff(x) * (head(y,-1) + tail(y,-1)))/2
	}
    }
    
    if(retx) {c(x,y)} else {round(y, 3)}
}


de.rocr = function(data, labels, measures='auc'){

    # Calculate AUC of ROC curve and average difference
    scores = sapply(rownames(data), function(a){calc_rocr(as.numeric(data[a,]), labels, measures=measures)})

    # Return marker genes
    markers = data.frame(gene=rownames(data), stringsAsFactors=FALSE)
    cname = paste(measures, collapse='_')
    markers[,cname] = scores
    markers = markers[order(markers[,cname], decreasing=T),]
    return(markers)
}


calc_fdr = function(predictions, labels){

    # Get classification cutoff with F measure
    f = calc_rocr(as.numeric(predictions), labels, 'f', retx=T)

    # Get true/false negative/positives
    u = factor(as.numeric(predictions >= f[[1]]), levels=c(0,1), ordered=T)
    q = as.vector(table(u, labels))
    tn = q[[1]]; fp = q[[2]]; fn = q[[3]]; tp = q[[4]];

    # Calculate statistics
    fdr = fp/(tp+fp)
    tpr = tp/(tp+fn)
    tnr = tn/(tn+fp)
    acc = (tp + tn)/sum(q)
    return(c(f[[1]], acc, tpr, tnr, fdr, f[[2]]))
}


de.fdr = function(data, labels, sens_cut=.1){

    # Calculate classification statistics
    markers = t(sapply(rownames(data), function(a){calc_fdr(data[a,], labels)}))
    colnames(markers) = c('cutoff', 'accuracy', 'sensitivity', 'specificity', 'fdr', 'f1')

    # Calculate average difference
    avg_diff = sapply(rownames(data), function(a){as.numeric(diff(tapply(as.numeric(data[a,]), labels, mean)))})

    # Return marker genes
    markers = cbind(markers, data.frame(avg_diff=avg_diff, gene=rownames(markers), stringsAsFactors=F))
    markers = markers[markers$sens >= sens_cut,]
    return(markers)

}


de.mast = function(data, covariates, formula=NULL, lrt_regex=TRUE, n.cores=1, glmer=FALSE){

    load_mast()
    options(mc.cores=n.cores)

    # Make single cell assay (SCA) object
    print(dim(data)); print(dim(covariates))
    fdata = data.frame(matrix(rep(1, nrow(data))))
    covariates = as.data.frame(covariates)
    sca = MAST::FromMatrix(as.matrix(data), covariates, fdata)

    # Fit MAST hurdle model
    if(is.null(formula)){
        formula = paste('~' , paste(colnames(covariates), collapse=' + '))
    }
    formula = as.formula(formula)
    if(glmer == FALSE){
        zlm.obj = zlm(formula, sca, force=TRUE)
    }
    
    # Likelihood ratio test
    if(is.logical(lrt_regex)){lrt_regex = colnames(covariates)}
    contrasts = grep(paste(lrt_regex, collapse='|'), colnames(zlm.obj@coefC), value=T, perl=T)
    print(contrasts)
    res = summary(zlm.obj, doLRT=contrasts)$datatable

    # Get component information
    res.f = res[res$component == 'logFC', .(primerid, contrast, coef)]
    res.d = res[res$component == 'D', .(primerid, contrast, coef, `Pr(>Chisq)`)]
    res.c = res[res$component == 'C', .(primerid, contrast, coef, `Pr(>Chisq)`)]
    res.h = res[res$component == 'H', .(primerid, contrast, `Pr(>Chisq)`)]

    # Combine results
    res = merge(res.d, res.c, by=c('primerid', 'contrast'), all=T, suffixes=c('D', 'C'))
    res = Reduce(function(...) merge(..., by=c('primerid', 'contrast'), all=T), list(res, res.f, res.h))
    res = data.frame(subset(res, !is.na(`Pr(>Chisq)`)), stringsAsFactors=F)

    # Cleanup results
    colnames(res) = c('gene', 'contrast', 'coefD', 'pvalD', 'coefC', 'pvalC', 'mastfc', 'pvalH')
    res = res[res$gene != '',]
    res = res[order(res$contrast, res$pvalH),]

    # Replace NAs in mastfc with maximum
    res = as.data.table(res)
    res[, mastfc := ifelse(is.na(mastfc), max(mastfc, na.rm=T), mastfc), .(contrast)]

    # Adjust p-values
    res[, padjD := p.adjust(pvalD, 'fdr'), .(contrast)]
    res[, padjC := p.adjust(pvalC, 'fdr'), .(contrast)]
    res[, padjH := p.adjust(pvalH, 'fdr'), .(contrast)]

    options(mc.cores=1)
    return(res)
}
