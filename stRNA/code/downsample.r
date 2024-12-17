
num_cells_per_group = function(groups, total_cells=NULL, cells_per_group=NULL){

    num_cells = sort(table(groups))

    if(!is.null(cells_per_group)){
        num_cells[num_cells > cells_per_group] = cells_per_group
    } else {
        n = sort(table(groups))
	if(length(n) == 1){
	    num_cells = total_cells
	    names(num_cells) = names(n)
	} else {
	    u = c(0, cumsum(n)[1:(length(n)-1)])
	    i = (total_cells - u)/seq(length(n), 1, -1) < n
	    if(sum(i) > 0){
	        num_cells[i] = as.integer(ceiling((total_cells - sum(n[!i]))/sum(i)))
	    }
	}
    }

    num_cells
}

resample = function(x,...){if(length(x)==1) x else sample(x,...)}

simple_downsample = function(cells, groups, ngene=NULL, total_cells=NULL, cells_per_group=NULL, replace=FALSE){

    # Downsample cells evenly across groups
    # Option: provide ngene to select HQ cells

    # Set ngene (default = 1)
    if(is.null(ngene)){
        ngene = structure(rep(1, length(cells)), names=cells)
    }
    if(is.null(names(ngene))){names(ngene) = cells}

    # Calculate group sizes
    groups = as.factor(groups)
    num_cells_per_group = num_cells_per_group(groups=groups, total_cells=total_cells, cells_per_group=cells_per_group)

    # Downsample cells within each group
    ds.cells = sapply(levels(groups), function(a){

        # Shuffle cells within group
        cells = resample(cells[groups == a], replace=replace)

	# Select by highest ngene
	cells[order(ngene[cells], decreasing=T)[1:num_cells_per_group[[a]]]]
    })
    ds.cells = as.character(na.omit(unname(unlist(ds.cells))))

    # Fix total cells
    if(!is.null(total_cells)){
	if(length(ds.cells) > total_cells){
            ds.cells = sort(resample(ds.cells, total_cells))
	}
    }

    return(ds.cells)
}

