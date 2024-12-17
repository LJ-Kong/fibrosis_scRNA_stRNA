select_bg_genes = function(obj=NULL, scores=NULL, genes.use, cells.use=NULL, nbins=20, nperm=100, type='mean', cutoff=0.25){

    # ---------------------------------------------------------------------------------------------------
    # Given singlecell object and target gene set, select [nperm] expression-matched background gene sets
    # ---------------------------------------------------------------------------------------------------
    # type:
    # - mean = match mean expression level across all cells
    # - num_ident = match number of idents with minimum expression cutoff (e.g. number of cell types with mean expression > 0.25)

    # fix input arguments
    if(is.null(cells.use)){
        cells.use = colnames(obj$data)
    }

    # divide genes into g equal-frequency expression bins
    if(is.null(scores)){
        if(type == 'mean'){
            scores = rowMeans(obj$data[,cells.use])
        } else {
            scores = nice_agg(t(obj$data), obj$ident)
	    scores = colSums(scores >= cutoff)
        }
    }
    genes2bin = setNames(cut2(scores, g=nbins), names(scores))
    bin2genes = sapply(levels(genes2bin), function(a) names(genes2bin)[genes2bin == a])

    # select control genes
    genes.use = intersect(genes.use, names(scores))

    # background gene sets
    genes.bg = sapply(genes.use, function(gene){
        bin = genes2bin[[gene]]
	sample(bin2genes[[bin]], nperm, replace=TRUE)
    })

    if(nperm > 1){genes.bg = t(genes.bg)}
    return(genes.bg)
}
