# Run fast PCA on a singlecell object
run_rpca = function(obj, data=NULL, k, genes.use=NULL, cells.use=NULL, rescale=FALSE, robust=FALSE, scale.max=Inf){

    # This function calculates a PCA on the DGE in obj$scale.data
    # If data, this function will still update the singlecell object
    # If cells.use, you may want to rescale data with rescale=TRUE

    # Get data and subset
    if(is.null(data)){data = obj$scale.data}
    if(is.null(genes.use)){genes.use = rownames(data)}
    if(is.null(cells.use)){cells.use = colnames(data)}
    data = data[genes.use, cells.use]

    # Calculate PCA
    if(robust == FALSE){
        print(paste('Calculating rpca on [cells x genes] matrix with center =', rescale, 'and scale =', rescale))
	data = t(data)
	if(rescale == TRUE){
	    data = scale(data)
	    data[data < -1*scale.max] = -1*scale.max
	    data[data > scale.max] = scale.max
	}
	pca.obj = rpca(data, center=FALSE, scale=FALSE, retx=TRUE, k=k)
    } else {
        print(dim(data))
        #pca.obj = rrpca(t(data), center=rescale, scale=rescale, retx=TRUE, k=k)
	pca.obj = rrpca(t(data), k=k)
    }
    obj$pca.obj = list(pca.obj)
    obj$pca.rot = data.frame(pca.obj$x)
    obj$pca.x = data.frame(pca.obj$rotation)

    # Return singlecell object
    return(obj)
}

