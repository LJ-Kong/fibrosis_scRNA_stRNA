
mem_cbind = function(M, cleanup=FALSE, verbose=FALSE){
    
    # Initialize matrix
    x = M[[1]]
    
    # Fast merge
    for(i in 2:length(M)){if(verbose){cat('.')}
        m = M[[i]]
	
	# Add rows to x
	r = setdiff(rownames(m), rownames(x))
	if(length(r) > 0){
	    n = Matrix(0, nrow=length(r), ncol=ncol(x), sparse=T)
	    rownames(n) = r
	    x = rbind(x, n)
	}
	
	# Add rows to m
	r = setdiff(rownames(x), rownames(m))
	if(length(r) > 0){
	    n = Matrix(0, nrow=length(r), ncol=ncol(m), sparse=T)
	    rownames(n) = r
	    m = rbind(m, n)
	}
	
	# Align rows
	m = m[rownames(x),,drop=F]
	
	# Combine data
	x = cbind(x, m)

	# Cleanup data
	if(cleanup){
	    M[[i]] = NA
	}
    }
    
    x
}

mem_rbind = function(M, cleanup=FALSE){
    
    # Initialize matrix
    x = M[[1]]
    
    # Fast merge
    for(i in 2:length(M)){
        m = M[[i]]
	
	# Add columns to x
	r = setdiff(colnames(m), colnames(x))
	if(length(r) > 0){
	    n = Matrix(0, ncol=length(r), nrow=nrow(x), sparse=T)
	    colnames(n) = r
	    x = cbind(x, n)
	}
	
	# Add columns to m
	r = setdiff(colnames(x), colnames(m))
	if(length(r) > 0){
	    n = Matrix(0, ncol=length(r), nrow=nrow(m), sparse=T)
	    colnames(n) = r
	    m = cbind(m, n)
	}
	
	# Align columns
	m = m[,colnames(x),drop=F]
	
	# Combine data
	x = rbind(x, m)
	
	# Cleanup data
	if(cleanup){
	    M[[i]] = NA
	}
    }
    
    x
}


sparse_cbind = function(M, verbose=FALSE){
    
    # Initialize matrix
    rows = unique(sort(unlist(sapply(M, rownames))))
    x = Matrix(0, nrow=length(rows), sparse=T)
    rownames(x) = rows
    
    # Fast merge with rbind/cbind
    for(m in M){if(verbose){cat('.')}
        
        # Align rows
	if(length(rows) - nrow(m) > 0){
	    n = Matrix(0, nrow=length(rows)-nrow(m), ncol=ncol(m), sparse=T)
	    rownames(n) = setdiff(rows, rownames(m))
	    m = rbind(m, n)
	}
	m = m[rows,,drop=F]
	
	# Combine data
	x = cbind(x, m)
    }

    # Remove first column
    x = x[,2:ncol(x),drop=F]
    x
}


sparse_rbind = function(M){
    
    # Initialize matrix
    cols = unique(sort(unlist(sapply(M, colnames))))
    x = Matrix(0, ncol=length(cols), sparse=T)
    colnames(x) = cols

    # Fast merge with rbind/cbind
    for(m in M){

        # Align columns
	n = Matrix(0, ncol=length(cols)-ncol(m), nrow=nrow(m), sparse=T)
	colnames(n) = setdiff(cols, colnames(m))
	m = cbind(m, n)
	m = m[,cols,drop=F]

	# Combine data
	x = rbind(x, m)
    }

    # Remove first row
    x = x[2:nrow(x),,drop=F]
    x
}

sparse_matsum = function(M){
    
    # Initialize matrix
    rows = unique(sort(unlist(sapply(M, rownames))))
    cols = unique(sort(unlist(sapply(M, colnames))))
    x = Matrix(0, nrow=length(rows), ncol=length(cols), sparse=T)
    rownames(x) = rows
    colnames(x) = cols

    # Fast merge with rbind/cbind
    for(m in M){        
	i = intersect(rownames(x), rownames(m))
	j = intersect(colnames(x), colnames(m))
	x[i,j,drop=F] = x[i,j,drop=F] + m[i,j,drop=F]
    }
    x
}


fast_rnorm = function(x){
    scaleMargins(x, rows=1/rowSums(x))
}


fast_cnorm = function(x){
    scaleMargins(x, cols=1/colSums(x))
}


mtx_filenames = function(prefix, data='matrix.mtx', rows='genes.tsv', cols='barcodes.tsv'){
    
    # Get sparse matrix filenames from prefix and patterns
    
    # Get prefix
    if(dir.exists(prefix)){prefix = paste0(prefix, '/')}
    if(!file.exists(prefix)){prefix = paste0(prefix, '.')}
    
    # Fix names
    prefix = gsub('//', '/', prefix)
    #prefix = gsub('\\.\\.', '.', prefix)
    prefix = gsub('/\\.', '/', prefix)
    
    # Get filenames
    if(!file.exists(data)){data = paste0(prefix, data)}
    if(!file.exists(rows)){rows = paste0(prefix, rows)}
    if(!file.exists(cols)){cols = paste0(prefix, cols)}
    
    # Check gunzip extension
    if(!file.exists(data)){data = paste0(data, '.gz')}
    if(!file.exists(rows)){rows = paste0(rows, '.gz')}
    if(!file.exists(cols)){cols = paste0(cols, '.gz')}

    # Remove gunzip extension
    if(!file.exists(data)){data = gsub('\\.gz', '', data)}
    if(!file.exists(rows)){rows = gsub('\\.gz', '', rows)}
    if(!file.exists(cols)){cols = gsub('\\.gz', '', cols)}

    return(list(data=data, rows=rows, cols=cols))
}


read_mtx = function(prefix, data='matrix.mtx', rows='features.tsv', cols='barcodes.tsv', filter=FALSE, fix_duplicates=FALSE){

    # Get filenames
    fns = mtx_filenames(prefix=prefix, data=data, rows=rows, cols=cols)
    print(fns)
    
    # Fix 10X inconsistencies
    if(!all(sapply(fns, file.exists))){
        fns = mtx_filenames(prefix=prefix, data=data, rows='genes.tsv', cols='barcodes.tsv')
    }
    data = fns$data
    rows = fns$rows
    cols = fns$cols
    
    # Read data
    print(paste('Reading', data))
    data = readMM(data)
    print(dim(data))
    rows = read.table(rows, sep='\t', stringsAsFactors=F, quote='', comment.char='')
    
    # Fix 10X inconsistencies
    if(ncol(rows) > 1){
        rows = rows[,sapply(rows, function(a) ! all(a == 'Gene Expression')),drop=F]
	rows = rows[,ncol(rows)]
    } else {
        rows = rows[,1]
    }
    print(length(rows))
    
    cols = read.table(cols, stringsAsFactors=F)
    if(ncol(cols) > 1){
        cols = cols[,ncol(cols)]
    } else {
        cols = cols[,1]
    }
    print(length(cols))
    
    # Filter matrix
    if(filter == TRUE){
        print('read_mtx: minimal filtering (requiring at least 10 genes per cell)')
	j = colSums(data > 0) >= 10
	data = data[,j]
	cols = cols[j]
    }
    
    # Set names
    if(fix_duplicates == TRUE){
        print('read_mtx: fixing duplicate rownames')
	data = as(rowsum(as.matrix(data), rows), 'sparseMatrix')
    } else {
        rownames(data) = rows
    }
    colnames(data) = cols
    
    return(data)
}


read_mtx2 = function(prefix, filter=FALSE, fix_duplicates=FALSE){
    
    # Get filenames
    fns = list.files(pattern=prefix)
    data = grep('matrix', fns, value=T)
    rows = grep('genes|features', fns, value=T)
    cols = grep('barcodes', fns, value=T)
    
    # Read data
    print(paste('Reading', data))
    data = readMM(data)
    if(grepl('csv', rows)){rsep=','} else {rsep='\t'}
    if(grepl('csv', cols)){csep=','} else {csep='\t'}
    rows = read.table(rows, sep=rsep, stringsAsFactors=F, quote='', comment.char='', header=F)
    cols = read.table(cols, sep=csep, stringsAsFactors=F, quote='', comment.char='', header=F)
    
    # Select column with most unique entries
    if(ncol(rows) > 1){
        rows = rows[,which.max(apply(rows, 2, function(a) length(unique(a))))]
    }
    if(ncol(cols) > 1){
        cols = cols[,which.max(apply(cols, 2, function(a) length(unique(a))))]
    }

    # Switch rows and columns if needed
    if(abs(length(rows) - nrow(data)) > (length(cols) - nrow(data))){
        print('Switching rows and columns')
	data = t(data)
    }
    
    # Smart alignment
    if((length(rows) - 1) == nrow(data)){
        rows = rows[2:length(rows)]
    }
    if((length(cols) - 1) == ncol(data)){
        cols = cols[2:length(cols)]
    }    
    
    # Filter matrix
    if(filter == TRUE){
        print('read_mtx: minimal filtering (requiring at least 10 genes per cell)')
	j = colSums(data > 0) >= 10
	data = data[,j]
	cols = cols[j]
    }
    
    # Set names
    if(fix_duplicates == TRUE){
        print('read_mtx: fixing duplicate rownames')
	data = as(rowsum(as.matrix(data), rows), 'sparseMatrix')
    } else {
        rownames(data) = rows
    }
    colnames(data) = cols
    
    return(data)
}


write_mtx = function(x, prefix='.', data='matrix.mtx', rows='genes.tsv', cols='barcodes.tsv', temp=FALSE){

    # Get filenames
    if(temp == TRUE){prefix = tempfile(tmpdir='~/tmp')}
    fns = mtx_filenames(prefix=prefix, data=data, rows=rows, cols=cols)
    data = fns$data
    rows = fns$rows
    cols = fns$cols
    
    # Write data
    writeMM(x, data)
    writeLines(rownames(x), rows)
    writeLines(colnames(x), cols)

    # Return filenames
    list(data=data, rows=rows, cols=cols)
}


colMeans_drop0 = function(x){
  n = diff(x@p)
  n[n == 0] = 1
  setNames(colSums(x) / n, 1:ncol(x))
}


rowMeans_drop0 = function(x){
  i = x@i + 1
  n = tabulate(i)
  n[n == 0] = 1
  setNames(rowSums(x) / n, 1:nrow(x))
}


align_mtx = function(x, rows.use=NULL, cols.use=NULL, by_row=TRUE, by_col=TRUE, missing=NA){

    # align a list of matrices by rows and/or columns
    # -----------------------------------------------

    # get row and column names
    if(is.null(rows.use)){
        rows.use = sort(unique(unlist(lapply(x, rownames))))
    }
    if(is.null(cols.use)){
        cols.use = sort(unique(unlist(lapply(x, colnames))))
    }
    
    # check if all matrices are sparse
    sparse = all(sapply(x, function(a) is(a, 'sparseMatrix')))
    
    # align matrices along rows/columns
    y = lapply(x, function(xi){
        if(! by_row){rows.use = rownames(xi)}
	if(! by_col){cols.use = colnames(xi)}
	if(sparse == TRUE){
	    yi = Matrix(missing, nrow=length(rows.use), ncol=length(cols.use), sparse=TRUE)
	    rownames(yi) = rows.use
	    colnames(yi) = cols.use
	    yi[rownames(xi), colnames(xi)] = as(xi, 'sparseMatrix')
	} else {
	    yi = matrix(missing, nrow=length(rows.use), ncol=length(cols.use))
	    rownames(yi) = rows.use
	    colnames(yi) = cols.use
	    yi[rownames(xi), colnames(xi)] = as.matrix(xi)
	}
	yi
    })
    y
}
