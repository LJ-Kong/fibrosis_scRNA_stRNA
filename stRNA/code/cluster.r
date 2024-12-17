cosine_dist = function(x){as.dist(1 - crossprod(x)/crossprod(t(matrix(sqrt(colSums(x**2))))))}

cross_nng = function(x, groups, k=10){
    groups = as.factor(groups)

    # Calculate "cross" nearest neighbor graph
    dx = lapply(levels(groups), function(g1){print(g1)
  	i = which(groups == g1)
        u = lapply(levels(groups), function(g2){
	    j = which(groups == g2)
	    dx = get.knnx(data=x[j,], query=x[i,], k=k)$nn.index
	    dx = apply(dx, 2, function(a) j[a])
	    dx
	})
	u = do.call(cbind, u)
        rownames(u) = rownames(x)[i]
	u
    })
    dx = do.call(rbind, dx)
    dx = dx[rownames(x),]

    # This code is from the "nng" function
    edges = matrix(unlist(sapply(1:nrow(x), function(i) {
        rbind(rep(i, k), dx[i, ])
    })), nrow = 2)
    n =  nrow(x)
    graph(edges, n = n, directed = TRUE)
}

run_graph_cluster = function(data, k=100, groups=NULL, method='infomap', weighted=FALSE, dist='cosine', do.fast=FALSE, out=NULL, knn=NULL, algo=c("kd_tree", "cover_tree", "CR", "brute")){

    if(is.null(knn)){

        print('Building kNN graph')
	if(is.null(groups)){
            knn = nng(data, k=k, method=dist, use.fnn=do.fast, algorithm=algo)
	} else {
	    new_k = as.integer(k/length(unique(groups)))
	    print(paste0('Calculating cross-kNN with ', length(unique(groups)), ' groups and k =', new_k))
	    knn = cross_nng(data, groups=groups, k=k)
	}
        if(weighted == TRUE){

	    print('Calculating Jaccard similarity')
            s = similarity(knn, method='jaccard')

	    print('Building weighted graph')
            knn = graph.adjacency(s, mode='undirected', weighted=T)

        }
    }

    if(method == 'louvain'){
        print('Louvain clustering')
        m = cluster_louvain(as.undirected(knn))$membership
    }

    if(method == 'infomap'){
        print('Infomap clustering')
        m = cluster_infomap(knn)$membership
    }

    # Write output
    print(sprintf('Clustering with k = %d finished', k))
    if(!is.null(out)){
        write.table(m, file=out, sep='\t', quote=F)
    }

    # Return clusters
    clusters = data.frame(x=as.character(m), row.names=rownames(data), stringsAsFactors=F)
    colnames(clusters) = c(paste0(method, '.', dist, '.k', k))
    return(clusters)
}

