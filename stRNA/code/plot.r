
qtrim = function(x, qmin=0, qmax=1, vmin=-Inf, vmax=Inf, rescale=NULL){

    # Trim by value
    x[x < vmin] = vmin
    x[x > vmax] = vmax

    # Trim by quantile
    u = quantile(x, qmin, na.rm=T)
    v = quantile(x, qmax, na.rm=T)
    x[x < u] = u
    x[x > v] = v

    return(x)
}

subsample_points = function(coords, k, nbins=25, bin_type='size'){

    # Evenly subsample points across multiple axes
    # - coords = data.frame(x1=axis1, x2=axis2, ...)
    # - k = number of points to subsample
    # - nbins = number of bins per axis
    # - bin_type = 'size' (equal size) or 'freq' (equal frequency)
    # For example, splits an xy plot into boxes and samples points from each box
    # Return: TRUE/FALSE vector for each row of coords

    # Divide points into bins
    if(bin_type == 'size'){
        g = interaction(lapply(coords, cut, breaks=nbins))
    } else if(bin_type == 'freq'){
        g = interaction(lapply(coords, cut2, g=nbins))
    } else {
        print("Error: ! bin_type %in% c('size', 'freq')")
    }

    # Subsample points from each bin
    i = as.integer(simple_downsample(1:nrow(coords), groups=g, total_cells=k))
    1:nrow(coords) %in% i
}

load_signature = function(file=NULL, file.regex=NULL, file.cols=NULL, sep=''){
    sig = read.table(file, stringsAsFactors=F, row.names=1, sep=sep, quote='')
    sig = structure(strsplit(sig[,1], ','), names=rownames(sig))
    if(!is.null(file.regex)){file.cols = grep(file.regex, names(sig), value=T)}
    if(!is.null(file.cols)){sig = sig[file.cols]}
    return(sig)
}

plot_images = function(obj, ..., image.use=NULL, titles.use=NULL, do.image=TRUE, ncol=NULL, image_ds=1, add_image=TRUE, vmin=NULL, vmax=NULL, legend_scale=1, legend_width=.1, do.legend=T){
    sink('/dev/null')
    on.exit(if(sink.number() > 0){sink()})
    if(is.null(titles.use)){
        titles.use = setNames(names(obj$image), names(obj$image))
    }
    if(is.null(image.use)){
        image.use=names(obj$image)
    }
    ps = sapply(image.use, function(a){
        if(do.image == FALSE){
	    image.use = NULL
	} else {
	    image.use = obj$image[[a]]
	}
	if(add_image == TRUE){
            add_image = obj$image[[a]]
	} else {
	    add_image = NULL
	}
	plot_tsne(obj, coords=obj$impos[[a]], theme_void=T, do.label=F, image=image.use, image_ds=image_ds, add_image=add_image, vmin=vmin, vmax=vmax, ...) + ggtitle(titles.use[[a]])
    }, simplify=F)
    sink()
    if(is.null(vmin) | is.null(vmax)){print('WARNING: shared axis, but vmin/vmax not set')}
    legend = get_legend(ps[[1]])
    p = plot_grid(plotlist=lapply(ps, function(a) a + theme(legend.position='none', axis.text.x=element_blank(), axis.text.y=element_blank()) + xlab('') + ylab('')), ncol=ncol, align='hv')
    if(do.legend){p = plot_grid(p, legend, rel_widths=c(1-legend_width, legend_width))}
    p
}

plot_tsne = function(obj=NULL, names=NULL, scores=NULL, coords=NULL, data=NULL, meta=NULL, regex=NULL, files=NULL, file.cols=NULL, file.regex=NULL, top=NULL, ident=TRUE, data.use='log2',
                     combine_genes='mean', cells.use=NULL, ymin=0, ymax=1, num_col='auto', pt.size=.75, font.size=11, do.label=T, label.size=5, do.title=TRUE, title.use=NULL, title.size=10, dpal=NULL, cpal=NULL,
	             do.legend=TRUE, legend.title='', share_legend=FALSE, share_legend_rescale=FALSE, legend_width=.05, legend.cols=NULL, vmin=NA, vmax=NA,
		     na.value='transparent', out=NULL, nrow=1.5, ncol=1.5, return_plotlist=FALSE, ds_cells=NULL, agg.ident=FALSE, do.sort=FALSE, do.revsort=FALSE, image=NULL, image_ds=10, xzoom=NULL, yzoom=NULL, imborder=50,
		     do.density=FALSE, xlab='Axis 1', ylab='Axis 2', theme_void=FALSE,  add_image=NULL, multichannel=NULL, multichannel_rescale=TRUE, s.use=NULL,
		     contour=NULL, contour.size=.75, contour.color='white', contour.bins=2,
		     legend.key.size=1.2, legend.text.size=8, legend.title.size=9, cblind=F, ...){


    # TSNE coordinates
    d = structure(obj$tsne.rot[,1:2], names=c('x', 'y'))
    
    # Cell identities
    if(!is.logical(ident)){
        d$Identity = ident
    } else if(ident & !is.null(obj)){
        d$Identity = obj$ident
    }
        
    # Use alternative coordinates
    if(!is.null(coords)){
        if(is.null(rownames(coords))){
	    rownames(coords) = colnames(obj$data)
	}
	i = intersect(rownames(d), rownames(coords))
	d = d[i,]
	d[rownames(coords),c('x','y')] = coords[,1:2]
    }
        
    # Fix zoom coordinates
    if(!is.null(xzoom)){
        xzoom[1] = max(1, xzoom[1])
	xzoom[2] = min(dim(image)[2], xzoom[2])
    }
    if(!is.null(yzoom)){
        yzoom[1] = max(1, yzoom[1])
	yzoom[2] = min(dim(image)[1], yzoom[2])
    }
        
    # Cell scores
    if(!is.null(cells.use)){cells.use = intersect(cells.use, rownames(d))}
    scores = score_cells(obj=obj, data=data, meta=meta, names=names, regex=regex, files=files, file.cols=file.cols, file.regex=file.regex, top=top, scores=scores,
                         data.use='log2', combine_genes=combine_genes, cells.use=cells.use)
    
    if(!is.null(scores)){
        i = intersect(rownames(d), rownames(scores))
	ni = c(names(d), names(scores))
	d = cbind.data.frame(d[i,], scores[i,], stringsAsFactors=F)
	colnames(d) = ni
    }
    
    if(agg.ident == TRUE){
        j = setdiff(colnames(d), c('x', 'y', 'Identity'))
        d[,j] = apply(d[,j,drop=F], 2, function(a){
	    ave(a, obj$ident, FUN=function(b){
	        if(is.character(b) | is.factor(b)){names(sort(table(b), dec=T))[[1]]} else{mean(b, na.rm=T)}
	    })
	})
    }

    # Subset cells
    if(is.null(cells.use)){
        cells.use = rownames(d)
    }
    if(!is.null(ds_cells)){
        print(paste('Downsampling', ds_cells, 'cells evenly across xy-grid'))
        i = subsample_points(coords=d[,1:2], k=ds_cells)
	cells.use = cells.use[i]
	print(paste('Selected', length(cells.use), 'cells'))
    }

    if(is.null(cells.use)){cells.use = rownames(d)}
    d = data.frame(d[cells.use,])

    # Initialize plotlist
    cat('\nPlotting:', paste(colnames(subset(d, select=-c(x,y))), collapse=', '), '\n')
    ps = list()

    # Shuffle point order
    d = d[sample(1:nrow(d)),]
    
    # Create column for plotting image only
    if(!is.null(add_image)){d = cbind("H&E"=1, d); image=add_image}
    
    # Get limits for shared legend
    if(share_legend_rescale == TRUE){
        j = names(which(sapply(subset(d, select=-c(x,y)), is.numeric)))
	cat('\nShared limits:', paste(j, collapse=', '), '\n')
        if(is.na(vmin)){vmin = na.omit(min(qtrim(d[,j], qmin=ymin, qmax=ymax, vmin=vmin, vmax=vmax)))}
	if(is.na(vmax)){vmax = na.omit(max(qtrim(d[,j], qmin=ymin, qmax=ymax, vmin=vmin, vmax=vmax)))}
	cat('> vmin =', vmin, '\n> vmax =', vmax, '\n')
    }
    
    if(!is.null(multichannel)){
        
        # convert multichannel data to rgb values
	m = as.data.frame(multichannel)
	if(is.null(image)){i = intersect(rownames(d), rownames(m)); d=d[i,]; m=m[i,]}
	if(ncol(m) == 2){m$Blank = 0}
	l = factor(colnames(m), levels=colnames(m))
	o = apply(m, 1, max)
	o[is.na(o)] = 0
	if(multichannel_rescale){print('rescaling values for multiscale plot')
	    m = scale(m, center=F, scale=apply(m, 2, quantile, ymax))
	}
	m[is.na(m)] = 0
	m[m > 1] = 1
	m = apply(m, 1, function(a) rgb(a[[1]], a[[2]], a[[3]], maxColorValue=1))

	# save rgb to "multichannel" column of d
	d$multichannel = NA
	i = intersect(rownames(d), names(m))
	d[i,'multichannel'] = m[i]
	col.order = c('x','y','multichannel', setdiff(colnames(d), c('x','y','multichannel')))
	d = d[,col.order]
	
	# re-order with do.sort
	o[setdiff(rownames(d), names(o))] = 0
	o = names(sort(o))
	d = d[o,]
    }

    if(!is.null(contour)){
	d$contour = NA
	i = intersect(rownames(d), names(contour))
	d[i,'contour'] = contour[i]
    }
    
    # easy titles
    col.use = setdiff(colnames(d), c('x', 'y', 'contour'))
    if(length(col.use) == length(title.use)){
        title.use = setNames(title.use, col.use)
    }
    
    # select image to use
    if(!is.null(s.use)){
        if(s.use == 'maxmean'){
	    s.use = names(which.max(sapply(vis$impos, function(a){
	        i = rownames(a)
		j = setdiff(colnames(d), c('x', 'y', 'Identity'))
		D = d[rownames(d) %in% i,j,drop=F]
		D = ifelse(is.na(D), NA, ifelse(is.numeric(D), D, 1))
		mean(apply(D, 1, prod), na.rm=T)
	    })))
	    coords = obj$impos[[s.use]]
        } else if(s.use == 'maxsum'){
	    s.use = names(which.max(sapply(vis$impos, function(a){
	        i = rownames(a)
		j = setdiff(colnames(d), c('x', 'y', 'Identity'))
		D = d[rownames(d) %in% i,j,drop=F]
		D = ifelse(is.na(D), NA, ifelse(is.numeric(D), D, 1))
		sum(apply(D, 1, prod), na.rm=T)
	    })))
	    coords = obj$impos[[s.use]]
	} else if(s.use %in% names(obj$impos)){
	    coords = obj$impos[[s.use]]
	} else {
	    stop('error: ! s.use %in% c("maxmean", "maxsum", names(obj$impos))')
	}
        d = d[rownames(d) %in% rownames(coords),]
    }
    
    for(col in col.use){

        # plot NAs first
        d = d[c(which(is.na(d[,col])), which(!is.na(d[,col]))),]
	
	# order points by value
	if(do.sort == TRUE){if(!is.null(multichannel)){print('warning: re-sorting multichannel plot (set do.sort=F for multichannel)')}
	    if(!is.numeric(d[,col])){
		do = order(d[,col])
	    } else {
	        do = order(abs(d[,col]))
	    }
	    d = d[do,]
	}
	if(do.revsort == TRUE){
	    d = d[rev(order(d[,col])),]
	}	
	
	if(col == 'multichannel'){

	    # set up basic multichannel plot
            p = ggplot(d)	    
	    if(!is.null(image)){
	        p = p + plot_image(image, xzoom=xzoom, yzoom=yzoom, downsample=image_ds)
	    }
	    
	    # build multichannel legend
	    D = d[1:3,]
	    D$col = l
	    
	    p = p + geom_point(data=D, aes(x=x, y=y, col=col), stroke=0) + scale_colour_manual('', values=set.colors[c(1,3,2)])
	    
	    # plot multichannel rgb values
	    p = p + geom_point(aes_string(x='x', y='y'), colour=d[,col], size=pt.size, ...) +
	        theme_cowplot(font_size=font.size) + xlab(xlab) + ylab(ylab)
	
	} else if(is.numeric(d[,col])){

	    # Continuous plot
	    
	    # Get colors
	    if(!is.null(cpal)){
	        if(length(cpal) == 1){cont.colors = brewer.pal(9, cpal)} else {cont.colors = cpal}
	    } else {
	        cont.colors = material.heat(50)
	    }
	    if(cblind){cont.colors = viridis(15)}
	    d[,col] = qtrim(d[,col], qmin=ymin, qmax=ymax, vmin=vmin, vmax=vmax)
	    
	    p = ggplot(d)
	        
	    if(!is.null(image)){
	        p = p + plot_image(image, xzoom=xzoom, yzoom=yzoom, downsample=image_ds)
	    }
	    
	    if(is.null(add_image)){
	        p = p +
	        geom_point(aes_string(x='x',y='y',colour=col), size=pt.size, ...) +
		scale_colour_gradientn(colours=cont.colors, guide=guide_colourbar(barwidth=.5, title=legend.title), na.value=na.value, limits=c(vmin, vmax)) + 
	        theme_cowplot(font_size=font.size) + xlab(xlab) + ylab(ylab)
	    }

	} else {print('discrete')

	    # Discrete plot

	    # Get colors
	    if(do.label == TRUE){disc.colors = set2.colors} else {disc.colors = set2.colors}
	    
	    if(!is.null(dpal)){
	        if(is.null(names(dpal))){
		    # use alternative color palette
		    disc.colors = dpal
		} else {
		    # if named palette, then map each name to groups of "similar" colors
		    groups = sapply(names(dpal), function(a) grep(a, levels(d[,col]), value=T))
		    groups = groups[lengths(groups) > 0]
		    dpal = setNames(sapply(1:length(groups), function(i) nice_colors(n=length(groups[[i]]), col=dpal[[i]], type='single')), names(groups))
		    dpal = unlist(lapply(names(dpal), function(a){b = intersect(levels(d[,col]), groups[[a]]); setNames(dpal[[a]], b)}))
		    dpal = dpal[as.character(levels(d[,col]))]
		    disc.colors = dpal
		}
	    }
	    if(cblind){if(length(unique(d[,col])) <= 7){disc.colors = c("#D55E00", "#0072B2", "#009E73", "#CC79A7", "#56B4E9", "#E69F00", "#000000")}}
	    
	    if(! is.factor(d[,col])){d[,col] = factor(d[,col], levels=naturalsort(unique(d[,col])))}
	    p = ggplot(d)
	    
	    if(!is.null(image)){
	        p = p + plot_image(image, xzoom=xzoom, yzoom=yzoom, downsample=image_ds)
	    }

	    if(is.null(add_image)){
	    p = p +
	        geom_point(aes_string(x='x',y='y',colour=col), size=pt.size, ...) +
		theme_cowplot(font_size=font.size) +
		xlab('Axis 1') + ylab('Axis 2') +
		scale_colour_manual(values=disc.colors, na.value=na.value, drop=F) +
		theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
		guides(colour=guide_legend(ncol=legend.cols, title=legend.title))
	    }
	    
	    if(do.label == T){
	        t = aggregate(d[,c('x', 'y')], list(d[,col]), median)
		colnames(t) = c('l', 'x', 'y')
		p = p + geom_text_repel(data=t, aes(x=x, y=y, label=l, lineheight=.8), point.padding=NA, size=label.size, family='Helvetica') + theme(legend.position='none')
	    }
	}

	if(!is.null(contour)){
	    m = with(d[,c('x','y','contour')], interp(x, y, contour, nx=25, ny=25))
	    m = as.data.frame(cbind(expand.grid(m$x, m$y), as.numeric(m$z)))
	    colnames(m) = c('x','y','z')
	    p = p + geom_contour(data=m, aes(x=x, y=y, z=z), size=contour.size, color=contour.color, bins=contour.bins)
	}
	
	if(do.title == TRUE){
	    if(is.null(title.use)){title = col} else if(col %in% names(title.use)){title = title.use[[col]]} else {title = title.use}
	    p = p + ggtitle(title)
	}
	
	if(legend.title == ''){
	    p = p + theme(legend.title=element_blank())
	}
	if(do.legend == FALSE){
	    p = p + theme(legend.position='none')
	} else {
	    p = p + theme(legend.key.size=unit(legend.key.size, 'line'), legend.text=element_text(size=legend.text.size), legend.title=element_text(size=legend.title.size))
	}
	

	if(!is.null(image)){
	    if(!is.null(xzoom)){
	        ixmin = xzoom[1]
		ixmax = xzoom[2]
	    } else {
		ixmin = min(d$x) - imborder
		ixmax = max(d$x) + imborder
	    }
	    if(!is.null(yzoom)){
	        iymin=yzoom[1]
		iymax=yzoom[2]
	    } else {
		iymin = min(d$y) - imborder
		iymax = max(d$y) + imborder
	    }

	    p = p + xlim(c(ixmin, ixmax)) + ylim(c(iymin, iymax))
	}
	
	if(theme_void){p = p + theme_void()}
	if(do.title == TRUE){
	    p = p + theme(plot.title=element_text(size=title.size))
	} else {
	    p = p + theme(plot.title=element_blank())
	}
	    
	if(!is.null(add_image)){image=NULL; add_image=NULL; p = p + theme(legend.position='none')}
	
	ps[[col]] = p
    }

    if(num_col == 'auto'){
        num_col = ceiling(sqrt(length(ps)))
    }
    
    ps = make_compact(plotlist=ps, num_col=num_col)
    if(length(ps) > 1){
        if(share_legend == TRUE){
            p = share_legend(ps, num_col=num_col, width=legend_width)
        } else {
            p = plot_grid(plotlist=ps, ncol=num_col, align='hv')
        }
    }
    
    if(is.null(out)){
        p
    } else {
        save_plot(file=out, p, nrow=nrow, ncol=ncol)
    }
    
    if(return_plotlist == TRUE){return(ps)} else {return(p)}
}

plot_image = function(image, xzoom=NULL, yzoom=NULL, downsample=10){
    if(!is.null(xzoom)){
        xmin = xzoom[1]
	xmax = xzoom[2]
    } else {
        xmin = 1
	xmax = dim(image)[2]
    }
    if(!is.null(yzoom)){
        ymin = yzoom[1]
	ymax = yzoom[2]
    } else {
        ymin = 1
	ymax = dim(image)[1]
    }
    xdown = max(1, as.integer(downsample*(xmax-xmin)/dim(image)[2]))
    ydown = max(1, as.integer(downsample*(ymax-ymax)/dim(image)[1]))
    i = rev(seq(from=ymin, to=ymax, by=xdown))
    j = seq(from=xmin, to=xmax, by=ydown)
    annotation_raster(image[i,j,], xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, interpolate=TRUE)
}

make_compact = function(plotlist, num_col, xlab='b', ylab='l', title=NULL){

    # Make a "plot_grid" plotlist compact by removing axes from interior points

    # Calculate plot indices
    ntot = length(plotlist)
    ncol = min(ntot, num_col)
    indices = list(
        't' = 1:ncol,
	'r' = which(1:ntot %% ncol == 0),
	'b' = rev(1:ntot)[1:ncol],
	'l' = which(1:ntot %% ncol == 1),
	'a' = 1:ntot
    )
    if(ntot == 1){indices[['l']] = 1}

    # Get directions
    xlab = strsplit(xlab, '')[[1]]
    ylab = strsplit(ylab, '')[[1]]
    if(!is.null(title)){title = strsplit(title, '')[[1]]}

    # Fix xlab
    i = setdiff(1:ntot, unlist(indices[xlab]))
    plotlist[i] = lapply(plotlist[i], function(p) p + theme(axis.text.x=element_blank(), axis.title.x=element_blank()))

    # Fix ylab
    i = setdiff(1:ntot, unlist(indices[ylab]))
    plotlist[i] = lapply(plotlist[i], function(p) p + theme(axis.text.y=element_blank(), axis.title.y=element_blank()))

    # Fix title
    if(!is.null(title)){
        i = setdiff(1:ntot, unlist(indices[title]))
        plotlist[i] = lapply(plotlist[i], function(p) p + theme(plot.title=element_blank()))
    }

    return(plotlist)

}

share_legend = function(plotlist, num_col, rel_widths=NULL, width=.1, align='hv'){

    # align = 'hv' makes plots identical sizes (but adds whitespace)

    # Get first legend in plotlist
    i = min(which(sapply(plotlist, function(p) 'guide-box' %in% ggplotGrob(p)$layout$name)))
    cat(paste('\nUsing shared legend:', names(plotlist)[i], '\n'))
    legend = get_legend(plotlist[[i]])

    # Remove all legends
    plotlist = lapply(plotlist, function(p) p + theme(legend.position='none'))

    # Make combined plot
    if(is.null(rel_widths)){rel_widths = rep(1, length(plotlist))}
    p = plot_grid(plotlist=plotlist, ncol=num_col, align=align, rel_widths=rel_widths)
    p = plot_grid(p, legend, ncol=2, rel_widths=c(1-width, width))

    return(p)
}

ggheatmap = function(data, Rowv='hclust', Colv='hclust', xlab='', ylab='', xsec=FALSE, ysec=FALSE, xstag=FALSE, xstag_space=.15, ystag=FALSE, ystag_space=.15, title='', legend.title='', xbreaks=NULL, ybreaks=NULL,
                     pal='nmf', do.legend=TRUE, font_size=7, groups=NULL, discrete=FALSE, dpal=set.colors, trans='identity',
                     out=NULL, nrow=1.25, ncol=1.25, qmin=0, qmax=1, vmin=-Inf, vmax=Inf, symm=FALSE, hclust_met='complete', border='#cccccc', replace_na=NA,
		     labRow=NULL, labCol=NULL, ret.order=FALSE, ret.legend=FALSE, pvals=NULL, max_pval=1, pval_border='black', pval_width=.25, na.value='#cccccc', do.label=F, lab.use=NULL, lab.min=NA, lab.offset=0){

    # Scale values
    if(is.null(rownames(data))){rownames(data) = 1:nrow(data)}
    if(is.null(colnames(data))){colnames(data) = 1:ncol(data)}
    data[is.na(data)] = replace_na
    if(! discrete){
        data = qtrim(data, qmin=qmin, qmax=qmax, vmin=vmin, vmax=vmax)
    }
    # Convert to long format
    x = as.data.frame(data) %>% rownames_to_column('row') %>% gather(col, value, -row)
    if(discrete){
        x$value = as.factor(x$value)
    } else {
        x$value = as.numeric(x$value)
    }
    
    # Add labels
    if(do.label){
        if(is.null(lab.use)){
	    x$label = x$value + lab.offset
            if(!is.na(lab.min)){x$label[x$label < lab.min] = ''}
	} else {
	    l = as.data.frame(lab.use) %>% rownames_to_column('row') %>% gather(col, label, -row)
	    x = merge(x, l, by=c('row', 'col'))
	}
    }
    
    # Facet by group
    if(!is.null(groups)){
        # groups is a named list mapping each x or y value to a group
	ri = sum(x$row %in% names(groups))
	ci = sum(x$col %in% names(groups))
	if(ri > ci){
	    x$group = as.factor(groups[x$row])
	} else {
	    x$group = as.factor(groups[x$col])
	}
    }
    
    # Merge data and p-values
    if(is.null(pvals)){
        x$pval = border
    } else {
        if(ncol(pvals) == 3){
	    colnames(pvals) = c('row', 'col', 'pval')
	    pvals$row = as.character(pvals$row)
	    pvals$col = as.character(pvals$col)
	} else {
	    pvals = as.data.frame(pvals) %>% rownames_to_column('row') %>% gather(col, pval, -row)
	    pvals$pval = as.numeric(pvals$pval)
	}
	if(length(intersect(x$row, pvals$row)) == 0){
	    colnames(pvals) = c('col', 'row', 'pval')
	}
	x = as.data.frame(merge(as.data.table(x), as.data.table(pvals), by=c('row', 'col'), all.x=TRUE))
	x$pval = ifelse(x$pval <= max_pval, pval_border, NA)
	x$pval[is.na(x$pval)] = border
	x = x[order(x$pval),]
    }
    
    # Order rows
    if(length(Rowv) > 1){rowv = Rowv; Rowv = 'Rowv'} else {rowv = rev(rownames(data))}
    if(length(Colv) > 1){colv = Colv; Colv = 'Colv'} else {colv = colnames(data)}
    if(nrow(data) <= 2){Rowv = 'none'}
    if(ncol(data) <= 2){Colv = 'none'}
    if(Rowv == 'hclust'){
        rowv = rev(rownames(data)[hclust(dist(na.replace(data, 0)), method=hclust_met)$order])
    }
    if(Colv == 'hclust'){
        colv = colnames(data)[hclust(dist(t(na.replace(data, 0))), method=hclust_met)$order]
    }
    if(Rowv == 'none'){
        rowv = rev(rownames(data))
    }
    if(Colv == 'none'){
        colv = colnames(data)
    }
    if(Rowv == 'min'){
        rowv = rev(rownames(data)[order(apply(data, 1, which.min))])
    }
    if(Rowv == 'max'){
    	i = match(colnames(data)[apply(data, 1, which.max)], colv)
	rowv = rev(rownames(data)[order(i)])
    }
    if(Colv == 'min'){
        i = match(rownames(data)[apply(data, 2, which.min)], rowv)
	colv = rev(colnames(data)[order(i)])
    }
    if(Colv == 'max'){
	i = match(rownames(data)[apply(data, 2, which.max)], rowv)
	colv = rev(colnames(data)[order(i)])
    }
    Rowv = rowv
    Colv = colv
    
    # Set order of row and column labels
    x$row = factor(x$row, levels=Rowv)
    x$col = factor(x$col, levels=Colv)
    
    # Get odd/even indices
    if(length(Rowv) > 1){
        r1 = seq(1, length(Rowv), by=2)
        r2 = seq(2, length(Rowv), by=2)
    } else {
        r1 = r2 = 1
    }
    if(length(Colv) > 1){
        c1 = seq(1, length(Colv), by=2)
        c2 = seq(2, length(Colv), by=2)
    } else {
        c1 = c2 = 1
    }
        
    # Get plot data
    if(length(pal)==1){if(pal == 'nmf'){pal = rev(colorRampPalette(nmf.colors)(101))[10:101]} else {pal = colorRampPalette(brewer.pal(9, pal))(101)}}

    # Plot significant boxes last
    x = x[rev(order(x$pval != 'black')),]
    
    if(discrete == TRUE){
        x$value = as.factor(x$value)
	pal = set.colors
    }
        
    # Plot with geom_tile
    p = ggplot(x) +
        geom_tile(aes(x=as.numeric(col), y=as.numeric(row), fill=value), color=x$pval, size=pval_width) +
	labs(x=xlab, y=ylab, title=title, fill=legend.title) +
	theme_cowplot(font_size=font_size) +
	theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5), axis.line=element_blank())
    
    if(!is.null(groups)){
        p = p + facet_grid(group ~ ., scales='free', space='free')
    }
    
    # add labels
    if(do.label == TRUE){
        p = p + geom_text(aes(x=as.numeric(col), y=as.numeric(row), label=label))
    }
    
    # Set scale
    if(discrete == TRUE){
        p = p + scale_fill_manual(values=dpal)
    } else {
    if(is.infinite(vmin)){vmin = min(x$value)}
    if(is.infinite(vmax)){vmax = max(x$value)}
    limits = c(vmin, vmax)
    if(symm == TRUE){
        values = c(min(x$value), 0, max(x$value))
	p = p + scale_fill_gradientn(colours=pal, values=scales::rescale(values), na.value=na.value, trans=trans)
    } else {
        p = p + scale_fill_gradientn(colours=pal, limits=limits, na.value=na.value, trans=trans)
    }
    }
    # Secondary x-axis
    if(is.logical(xsec)){
        if(xsec == FALSE){
            p = p + scale_x_continuous(breaks=1:length(Colv), labels=Colv, expand=c(0,0))
        } else {
	    p = p + scale_x_continuous(breaks=c1, labels=Colv[c1], sec.axis=dup_axis(breaks=c2, labels=Colv[c2]), expand=c(0,0)) + theme(axis.text.x.top= element_text(angle=90, hjust=0, vjust=.5))
        }
    } else{
        # xsec = list(Colv -> new label)
        p = p + scale_x_continuous(breaks=1:length(Colv), labels=Colv, sec.axis=dup_axis(breaks=1:length(Colv), labels=xsec[Colv]), expand=c(0,0)) + theme(axis.text.x.top= element_text(angle=90, hjust=0, vjust=.5))
    }

    # Secondary y-axis
    if(is.logical(ysec)){
        if(ysec == FALSE){
            p = p + scale_y_continuous(breaks=1:length(Rowv), labels=Rowv, expand=c(0,0))
        } else {
	    p = p + scale_y_continuous(breaks=r1, labels=Rowv[r1], sec.axis=dup_axis(breaks=r2, labels=Rowv[r2]), expand=c(0,0))
        }
    } else {
        # xsec = list(Colv -> new label)
        p = p + scale_y_continuous(breaks=1:length(Rowv), labels=Rowv, sec.axis=dup_axis(breaks=1:length(Rowv), labels=xsec[Rowv]), expand=c(0,0)) + theme(axis.text.x.top= element_text(angle=90, hjust=0, vjust=.5))
    }

    if(!is.null(xbreaks)){
        p = p + scale_x_discrete(breaks=xbreaks)
    }
    if(!is.null(ybreaks)){
        p = p + scale_y_discrete(breaks=ybreaks)
    }
    
    # Add axes with ggrepel (use both sides to fit as many labels as possible)
    if(!is.null(labRow)){
        p = p + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
	if(!is.logical(labRow)){
	    lab.use = intersect(Rowv, labRow)
	    lab.l = ifelse(Rowv %in% lab.use[seq(from=1, to=length(lab.use), by=2)], Rowv, '')
	    lab.r = ifelse(Rowv %in% lab.use[seq(from=2, to=length(lab.use), by=2)], Rowv, '')
	    legend = get_legend(p)
	    p = p + theme(legend.position='none', plot.margin=margin(l=-10, unit='pt'))
	    axis.l = add_ggrepel_axis(plot=p, lab=lab.l, type='left', font_size=font_size)
	    axis.r = add_ggrepel_axis(plot=p, lab=lab.r, type='right', font_size=font_size)
	    p = plot_grid(axis.l, p, axis.r, nrow=1, rel_widths=c(.15,1,.15), align='h', axis='tb')
	    p = plot_grid(p, legend, nrow=1, rel_widths=c(.925, .075))
	}
    }

    if(!is.null(labCol)){
        p = p + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
	if(xsec == TRUE & !is.logical(labCol)){
	    lab.use = intersect(Colv, labCol)
	    lab.u = ifelse(Colv %in% lab.use[seq(from=1, to=length(lab.use), by=2)], Colv, '')
	    lab.d = ifelse(Colv %in% lab.use[seq(from=2, to=length(lab.use), by=2)], Colv, '')
	    legend = get_legend(p)
	    p = p + theme(legend.position='none', plot.margin=margin(t=-10, b=-12.5, unit='pt'))
	    axis.u = add_ggrepel_axis(plot=p, lab=lab.u, type='up', font_size=font_size, do.combine=FALSE)
	    axis.d = add_ggrepel_axis(plot=p, lab=lab.d, type='down', font_size=font_size, do.combine=FALSE)
	    p = plot_grid(axis.u, p, axis.d, nrow=3, rel_heights=c(.15,1,.15), align='v', axis='lr')
	    p = plot_grid(p, legend, nrow=1, rel_widths=c(.925, .075))
	}
	if(xsec == FALSE & !is.logical(labCol)){
	    lab.use = ifelse(Colv %in% intersect(Colv, labCol), Colv, '')
	    legend = get_legend(p)
	    p = p + theme(legend.position='none')
	    axis.d = add_ggrepel_axis(plot=p, lab=lab.use, type='down', font_size=font_size, do.combine=FALSE)
	    p = plot_grid(p, axis.d, nrow=2, rel_heights=c(1, .15), align='v', axis='lr')
	    p = plot_grid(p, legend, nrow=1, rel_widths=c(.925, .075))
	}
    }
    
    if(xstag == TRUE){
        nudge = rep(c(0, 1), length.out=length(Colv))
	p = add_ggrepel_axis(plot=p, lab=Colv, dir='x', type='down', font_size=font_size, force=0, axis.width=xstag_space, nudge=nudge, ret.legend=ret.legend)
    }

    if(ystag == TRUE){
        nudge = rep(c(0, 1), length.out=length(Rowv))
	p = add_ggrepel_axis(plot=p, lab=Rowv, dir='y', type='left', font_size=font_size, force=0, axis.width=ystag_space, nudge=nudge, ret.legend=ret.legend)
    }

    if(do.legend == FALSE){p = p + theme(legend.position='none')}

    # Save plot
    if(!is.null(out)){
        save_plot(p, file=out, nrow=nrow, ncol=ncol)
    }

    if(ret.order == FALSE){p} else {list(p=p, Rowv=Rowv, Colv=Colv)}
}

add_ggrepel_axis = function(plot, lab, dir='both', type='left', lab.pos=NULL, lab.lim=NULL, nudge=0, force=1, axis.width=.10, legend.width=.075, font_size=8, ret.legend=FALSE, do.combine=TRUE){

    # Fix input arguments
    if(is.null(lab.pos)){lab.pos = 1:length(lab)}
    if(is.null(lab.lim)){lab.lim = c(min(lab.pos)-1, max(lab.pos)+1)}

    # Set label directions
    if(type == 'left'){x=0; y=lab.pos; xlim=c(-1,0); ylim=lab.lim; angle=0; nudge_x=-1*nudge; nudge_y=0}
    if(type == 'up'){x=lab.pos; y=0; xlim=lab.lim; ylim=c(0,1); angle=90; nudge_x=0; nudge_y=nudge}
    if(type == 'right'){x=0; y=lab.pos; xlim=c(0,1); ylim=lab.lim; angle=0; nudge_x=nudge; nudge_y=0}
    if(type == 'down'){x=lab.pos; y=0; xlim=lab.lim; ylim=c(-1,0); angle=90; nudge_x=0; nudge_y=-1*nudge}

    # Get data for ggplot
    d = data.frame(x=x, y=y, lab=lab, dir=dir, nudge_y=nudge_y)
    #d = d[((d$x %in% c(min(d$x), max(d$x))) | (d$y %in% c(min(d$y), max(d$y)))) | (d$lab != '')]
    d = d[d$lab != '',,drop=F]
        
    # Make ggrepel axis
    axis = ggplot(d, aes(x=x, y=y, label=lab)) +
           geom_text_repel(min.segment.length=grid::unit(0,'pt'), color='grey30', size=(font_size-1)/.pt, angle=angle, segment.color='#cccccc', segment.size=.15,
	                   direction=dir, nudge_x=nudge_x, nudge_y=nudge_y, force=force) +
	   scale_x_continuous(limits=xlim, expand=c(0,0), breaks=NULL, labels=NULL, name=NULL) +
	   scale_y_continuous(limits=ylim, expand=c(0,0), breaks=NULL, labels=NULL, name=NULL) +
    	   theme(panel.background = element_blank(), plot.margin = margin(0, 0, 0, 0, 'pt'))

    if(do.combine == FALSE){return(axis)}

    # Get plot legend
    legend = get_legend(plot)
    plot = plot + theme(legend.position='none')

    # Combine plots
    if(type == 'left'){
        plot = plot + scale_y_continuous(limits=lab.lim, expand=c(0,0))
        plot = plot + theme(plot.margin=margin(l=-12.5, unit='pt')) + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
	p = plot_grid(axis, plot, nrow=1, align='h', axis='tb', rel_widths=c(axis.width, 1-axis.width))
	if(ret.legend == FALSE){p = plot_grid(legend, p, nrow=1, rel_widths=c(legend.width, 1-legend.width))}
    }
    if(type == 'up'){
        plot = plot + scale_x_continuous(limits=lab.lim, expand=c(0,0))
        plot = plot + theme(plot.margin=margin(t=-10, unit='pt')) + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
	p = plot_grid(axis, plot, nrow=2, align='v', axis='lr', rel_heights=c(axis.width, 1-axis.width))
	if(ret.legend == FALSE){p = plot_grid(legend, p, nrow=1, rel_widths=c(legend.width, 1-legend.width))}
    }
    if(type == 'right'){
        plot = plot + scale_y_continuous(limits=lab.lim, expand=c(0,0))
        plot = plot + theme(plot.margin=margin(r=-10, unit='pt')) + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
	p = plot_grid(plot, axis, nrow=1, align='h', axis='tb', rel_widths=c(1-axis.width, axis.width))
	if(ret.legend == FALSE){p = plot_grid(p, legend, nrow=1, rel_widths=c(1-legend.width, legend.width))}
    }
    if(type == 'down'){
        plot = plot + scale_x_continuous(limits=lab.lim, expand=c(0,0))
        plot = plot + theme(plot.margin=margin(b=-11.5, t=5, unit='pt')) + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
	p = plot_grid(plot, axis, nrow=2, align='v', axis='lr', rel_heights=c(1-axis.width, axis.width))
	if(ret.legend == FALSE){p = plot_grid(p, legend, nrow=1, rel_widths=c(1-legend.width, legend.width))}
    }

    if(ret.legend == FALSE){p} else {list(p=p, legend=legend)}
}

simple_scatter = function(x, y, lab=NA, sig=FALSE, col=NULL, col.title='', size=NULL, size.title='', lab.use=NULL, lab.sig=FALSE, lab.near=0,  lab.n=0, lab.g=0, groups=NULL, lab.size=4, lab.type='up',
                          palette=NULL, dpal=NULL, cpal=NULL, legend.fontsize=10, border=F, edges=NULL, na.value='#cccccc', do.label=F,
                          xlab=NULL, ylab=NULL, out=NULL, nrow=1, ncol=1, min_size=1, max_size=3, xskip=c(0,0), yskip=c(0,0), xlim=NULL, ylim=NULL, alpha=1, unlab.grey=FALSE, auto_color='multi',
			  xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, xmin.n=1, xmax.n=1, ymin.n=1, ymax.n=1, cmin=-Inf, cmax=Inf,
			  xlog=FALSE, ylog=FALSE, xylog=FALSE, jitter=F, lm=FALSE, lmcol='#cccccc', line=F,
			  legend.names=c('', 'lab.n', 'lab.use', 'lab.near'), do.sort=FALSE, do.revsort=FALSE, max.overlaps=10){

    x = as.numeric(x)
    y = as.numeric(y)
    if(is.null(col)){col = rep('', length(x))}
    if(is.null(sig)){sig = rep(NA, length(x))}
    if(is.null(size)){size = rep(1, length(x))}
    if(is.character(col)){col = as.factor(col)}
    if(is.null(dpal)){dpal = palette}
    if(is.null(cpal)){cpal = palette}
    if(!is.null(dpal) | !is.null(cpal)){palette = TRUE}
    
    if(!is.null(edges)){colnames(edges) = c('x', 'y', 'xend', 'yend')}
    
    # adjust range
    xmin = max(na.omit(sort(x))[[xmin.n]], xmin)
    xmax = min(na.omit(sort(x, decreasing=T))[[xmax.n]], xmax)
    ymin = max(na.omit(sort(y))[[ymin.n]], ymin)
    ymax = min(na.omit(sort(y, decreasing=T))[[ymax.n]], ymax)
    x[x < xmin] = xmin
    x[x > xmax] = xmax
    y[y < ymin] = ymin
    y[y > ymax] = ymax

    if(!is.factor(col)){
        col = pmax(cmin, col)
        col = pmin(cmax, col)
    }
    
    d = data.frame(x=x, y=y, lab=lab, col=col, size=size, flag='', sig=sig, stringsAsFactors=FALSE)
    i.lab = !is.na(d$lab)
    di = d[i.lab,]

    if(do.sort){d = d[order(d$col),]}
    if(do.revsort){d = d[order(-d$col),]}

    if(is.null(xlab)){xlab = deparse(substitute(x)); if((length(xlab) > 5) || (substr(xlab, 1, 2) == 'c(')){xlab=''}}
    if(is.null(ylab)){ylab = deparse(substitute(y)); if((length(ylab) > 5) || (substr(ylab, 1, 2) == 'c(')){ylab=''}}

    if(lab.n > 0 | lab.g > 0){

        i = c()

	if('up' %in% lab.type | 'down' %in% lab.type){

            # get breaks
	    if(!is.null(groups)){groups.use = groups} else {
	        breaks = seq(from=min(di$x, na.rm=T), to=max(di$x, na.rm=T), length.out=min(20, lab.n))
	        groups.use = cut(di$x, breaks=breaks, include.lowest=TRUE)
		groups.use[xskip[[1]] < di$x & di$x < xskip[[2]]] = NA
	    }

	    # get cells
	    if('up' %in% lab.type){
	        ngene = ifelse(yskip[[1]] < di$y & di$y < yskip[[2]], NA, di$y)
	        i = c(i, as.numeric(simple_downsample(cells=1:nrow(di), groups=groups.use, ngene=ngene, total_cells=lab.n)))
		i = c(i, order(-1*di$y)[1:lab.g])
	    }
	    if('down' %in% lab.type){
	    	ngene = ifelse(yskip[[1]] < di$y & di$y < yskip[[2]], NA, -1*di$y)
	        i = c(i, as.numeric(simple_downsample(cells=1:nrow(di), groups=groups.use, ngene=ngene, total_cells=lab.n)))
		i = c(i, order(di$y)[1:lab.g])
	    }
	}

	if('right' %in% lab.type | 'left' %in% lab.type){

	    # get breaks
	    if(!is.null(groups)){groups.use = groups} else {
	        breaks = seq(from=min(di$y, na.rm=T), to=max(di$y, na.rm=T), length.out=min(20, lab.n))
	        groups.use = cut(di$y, breaks=breaks, include.lowest=TRUE)
		groups.use[yskip[[1]] < di$y & di$y < yskip[[2]]] = NA
	    }

	    # get cells
	    if('right' %in% lab.type){
	        i = c(i, as.numeric(simple_downsample(cells=1:nrow(di), groups=groups.use, ngene=di$x, total_cells=lab.n)))
		i = c(i, order(-1*di$x)[1:lab.g])
	    }

	    if('left' %in% lab.type){
	        i = c(i, as.numeric(simple_downsample(cells=1:nrow(di), groups=groups.use, ngene=-1*di$x, total_cells=lab.n)))
		i = c(i, order(di$x)[1:lab.g])
	    }
	}

	if('random' %in% lab.type){

	    i = c(i, sample(1:nrow(di), lab.n))
	}

	di[unique(i), 'flag'] = 'lab.n'
    }
    d[i.lab,] = di

    if(!is.null(lab.use)){

        # label neighbors
	if(lab.near > 0){
	    u = as.matrix(d[, c('x','y')])
	    v = as.matrix(d[lab %in% lab.use, c('x','y')])
	    i = unique(sort(unlist(apply(dist.matrix(u,v,skip.missing=TRUE,method='euclidean'), 2, order)[1:(lab.near + 1),])))
	    d$flag[i] = 'lab.near'
	}

        # label points
        d$flag[d$lab %in% lab.use] = 'lab.use'
    }

    if(lab.sig == TRUE){d$flag[d$sig == TRUE] = 'lab.use'; d$sig = FALSE}

    d$lab[d$flag == ''] = ''
    d = d[order(d$flag),]
    d$flag = factor(d$flag, levels=c('', 'lab.n', 'lab.use', 'lab.near'), ordered=T)

    if(auto_color == 'none'){d$col = ''}
    if(auto_color == 'bw' & all(col == '')){d$col = ifelse(d$flag == '', '', 'lab.n')}
    if(auto_color == 'multi' & all(col == '')){d$col = d$flag}

    # plot labeled points last
    d = d[order(d$flag),]
    i = is.na(d$col) | d$col == ''
    d = rbind(d[i,], d[!i,])

    if(unlab.grey == TRUE){d[d$lab == '',]$col = ''}

    # hack: fix legend names
    if(!is.numeric(d$col)){
        d$col = as.factor(d$col)
    	#levels(d$col) = legend.names
    }
    
    if(! jitter){position = 'identity'} else {position='jitter'}
    
    p = ggplot(d)

    if(!is.null(edges)){p = p + geom_segment(data=edges, aes(x=x,y=y,xend=xend,yend=yend), color='black')}
    
    if(border == FALSE){
        if(length(unique(size)) == 1){
	    p = p + geom_point(aes(x=x, y=y, col=col), size=d$size, alpha=alpha, position=position)
	} else {
            p = p + geom_point(aes(x=x, y=y, col=col, size=size), alpha=alpha, position=position)
	}
    } else {
        if(length(unique(size)) == 1){
	    p = p + geom_point(aes(x=x, y=y, fill=col, stroke=ifelse(sig == TRUE, stroke, 0)), size=size, alpha=alpha, pch=21, colour=border, position=position)
	} else {
	    p = p + geom_point(aes(x=x, y=y, fill=col, size=size, stroke=ifelse(sig == TRUE, stroke, 0)), alpha=alpha, pch=21, colour=border, position=position)
	}
    }
    p = p + geom_text_repel(data=d[(d$lab != '') | (1:nrow(d) %in% sample(1:nrow(d), min(nrow(d), 1000))),], aes(x=x, y=y, label=lab), size=lab.size, segment.color='grey', max.overlaps=max.overlaps) +
	theme_cowplot() + xlab(xlab) + ylab(ylab) + theme(legend.title=element_text(size=legend.fontsize), legend.text=element_text(size=legend.fontsize))
	
    if(all(size == 1)){p = p + scale_size(guide = 'none', range=c(min_size, max_size))} else {p = p + scale_size(name=size.title, range=c(min_size, max_size))}

    if(!is.null(xlim)){p = p + xlim(xlim)}
    if(!is.null(ylim)){p = p + ylim(ylim)}
    if(xlog | xylog){p = p + scale_x_log10()}
    if(ylog | xylog){p = p + scale_y_log10()}
    
    # label colors on plot
    if(do.label){
        t = aggregate(d[,c('x', 'y')], list(d[,'col']), median)
        colnames(t) = c('l', 'x', 'y')
        p = p + geom_text_repel(data=t, aes(x=x, y=y, label=l, lineheight=.8), point.padding=NA, size=lab.size, family='Helvetica') + theme(legend.position='none')
    }
        
    if(border != FALSE){
        if(!is.null(palette)){
            if(!is.numeric(d$col)){
	        p = p + scale_fill_manual(name=col.title, values=dpal, na.value=na.value, breaks=levels(as.factor(d$col)))
	    } else {
	        p = p + scale_fill_gradientn(name=col.title, colors=cpal, na.value=na.value)
	    }
        } else {
            if(!is.numeric(d$col)){
	        if(length(unique(d$col)) == 1){
		    v.use = '#000000'; p = p + theme(legend.position='none')
		} else if(length(unique(d$col)) == 2){
		    v.use = brewer.pal(3, 'Paired')
		} else {
		    v.use = set.colors
		}
	        p = p + scale_fill_manual(name=col.title, values=v.use) #+ theme(legend.position='none')
	    } else {
	        p = p + scale_fill_gradientn(name=col.title, colors=material.heat(100), na.value=na.value)
	    }
        }
    } else {
        if(!is.null(palette)){
	    if(!is.numeric(d$col)){
                p = p + scale_color_manual(name=col.title, values=dpal, na.value=na.value, breaks=levels(as.factor(d$col)))
	    } else {
	        p = p + scale_color_gradientn(name=col.title, colors=cpal, na.value=na.value)
	    }
        } else {
            if(!is.numeric(d$col)){
	    	if(length(unique(d$col)) == 1){
		    v.use = '#000000'; p = p + theme(legend.position='none')
		} else if(length(unique(d$col)) == 2){
		    v.use = brewer.pal(3, 'Paired')
		} else {
		    v.use = set.colors
		}
	        p = p + scale_color_manual(name=col.title, values=v.use) #+ theme(legend.position='none')
	    } else {
	        p = p + scale_color_gradientn(name=col.title, colors=material.heat(100), na.value=na.value)
	    }
        }
    }

    if(line){p = p + geom_line(aes(x=x, y=y), col=lmcol, size=.5, linetype='dashed')}
    if(lm){p = p + geom_smooth(aes(x=x, y=y), method='lm', se=F, col=lmcol, size=.5, linetype='dashed')}

    if(!is.null(out)){save_plot(plot=p, filename=out, nrow=nrow, ncol=ncol)}

    return(p)
}

