ffread = function(a, sep='\t', row.names=NULL, header=TRUE, as.dt=FALSE){
    if(grepl('gz', a)){
        h = strsplit(readLines(a, n=1), sep)[[1]]
        x = fread(paste0('zcat ', a), sep=sep, skip=1, header=F)
    } else {
        h = strsplit(readLines(a, n=1), sep)[[1]]
        x = fread(a, skip=1, sep=sep, header=F)
    }
    if(length(h) == (ncol(x)-1)){
        print('Setting header=T and row.names = 1')
        header = TRUE
        row.names = 1
    }
    print(paste('ffread: reading', a, 'with header =', header, 'and row.names=', row.names))
    if(nrow(x) == 1){
        if(row.names == 1){
            x = read.table(a, sep=sep, header=F, row.names=row.names, skip=1)
        }
    } else {
        x = data.frame(x, row.names=row.names)
    }
    if(length(h) == (ncol(x) + 1)){
        h = h[2:length(h)]
    }
    if(header == FALSE){
        if(length(h) == ncol(x)){
            x = rbind(h, x)
        } else {
            stop('error: file dimensions')
        }
    } else {
        if(length(h) == ncol(x)){
            colnames(x) = h
        } else if(length(h) == ncol(x)+1){
            colnames(x) = h[2:length(h)]
        } else {
            stop('error: file dimensions')
        }
    }
    if(as.dt == TRUE){x = as.data.table(x)}
    return(x)
}

