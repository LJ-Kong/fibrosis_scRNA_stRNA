
nice_network = function(graph=NULL, edges=NULL, symm=TRUE, self=FALSE,
                        layout='fruchtermanreingold', layout_weights=TRUE,
			size_legend='nodes', size_trans='identity', # 'log2', etc
                        node_colors=NULL, node_sizes=NULL, node_scale=c(5,5), node_limits=NULL, node_alpha=1, node_pal=NULL, node_padding=.5, node_strokes=0.25,
			edge_colors=NULL, edge_colors_order=NULL, edge_alpha=1, edge_scale=c(.5,.5), edge_limits=NULL, edge_pal=NULL, color_sign=FALSE, edge_padding=.5, linewidth_pal=1:6,
			node_guide='legend', edge_guide='legend', size_guide='legend', linetype_guide='legend', do.legend=TRUE,
			node_title='', edge_title='', size_title='', linetype_title='',
			lab.use=NULL, node_lab.size=4, edge_lab.size=3, segment.color='grey',
			ggtitle='', drop_levels=FALSE, out=NULL, nrow=1, ncol=1, ret.data=FALSE,
			do.arrow=FALSE, arrow.gap=.015, min_degree=NA, max.overlaps=10
			) {
    
    # plot network from graph or edgelist
    # -----------------------------------
    # graph = [m x n] matrix
    # edges = data.frame(1=source, 2=target, 3=weight, 4=color, 5=size, 6=label, 7=linetype, ...)
    # lab.use = list of node/edge labels to use
    # node_sizes = list(label=size)
    # node_colors = list(label=color)
    # edge_colors = [m x n] matrix
    # node_scale = c(min, max)  * scales specify the minimum and maximium size of the plotting symbol
    # edge_scale = c(min, max) 
    # node_limits = c(min, max) * limits specify the minimum and maximum size of the value being plotted
    # edge_limits = c(min, max)
    
    # set arrow parameters
    # --------------------
    if(do.arrow == FALSE){
        arrow = NULL
	arrow.gap = 0
    } else {
        arrow = arrow(length=unit(6, "pt"), type="closed")
	arrow.gap = arrow.gap
    }
    
    # construct edgelist
    # ------------------
    if(!is.null(graph)){
        graph = as.data.frame(as.matrix.data.frame(as.data.frame(graph)))
	edges = gather(graph %>% rownames_to_column('source'), target, weight, -source)
    }
    if(! all(c('source', 'target', 'weight') %in% colnames(edges))){
        stop('check format')
    }
    edges$weight = as.numeric(edges$weight)
    
    # fix edgelist
    # ------------
    edges = as.data.table(edges)
    edges$source = as.character(edges$source)
    edges$target = as.character(edges$target)
    
    if(symm == TRUE){
        i = apply(edges[,.(source, target)], 1, function(a) paste(a, collapse=' '))
	edges = edges[!duplicated(i)]
    }
    if(self == FALSE){
        edges = edges[source != target]
    }
    edges = edges[abs(weight) > 0]
    
    # node attributes
    # ---------------
    node_labels = sort(unique(c(edges$source, edges$target)))
    
    node_sizes = unlist(node_sizes) # recently added
    if(is.null(node_sizes)){
        node_sizes = setNames(rep(1, length(node_labels)), node_labels)
    }
    if(is.null(names(node_sizes))){
        names(node_sizes) = node_labels
    }
    
    node_colors = unlist(node_colors) # recently added
    if(is.null(node_colors)){
        node_colors = setNames(rep(1, length(node_labels)), node_labels)
    }
    if(is.null(names(node_colors))){
        names(node_colors) = node_labels
    }
    if(is.null(levels(node_colors))){
        node_colors = as.factor(node_colors)
    }
    node_strokes = unlist(node_strokes) # recently added
    if(length(node_strokes) == 1){
        node_strokes = setNames(rep(node_strokes, length(node_labels)), node_labels)
    }
    
    
    # edge attributes
    # ---------------
    if(!is.null(edge_colors)){
        edge_colors = as.data.frame(as.matrix.data.frame(as.data.frame(edge_colors)))
	edge_colors = gather(edge_colors %>% rownames_to_column('source'), target, color, -source)
	edges = merge(edges, edge_colors, by=c('source', 'target'))
    }
    if(! is.null(edge_colors_order)){
        edges$color = factor(edges$color, levels=edge_colors_order)
    }
    if(! 'color' %in% colnames(edges)){
        edges$color = 1
    }
    if(! 'label' %in% colnames(edges)){
        edges$label = ''
    }
    if(! 'linetype' %in% colnames(edges)){
        edges$linetype = 1
    }
    edges$color = as.factor(edges$color)
    edges$linetype = as.factor(edges$linetype)
    
    # fix edge weights
    # ----------------
    if(any(edges$weight < 0)){
        print('Warning: negative edges')
    }
    if(color_sign == TRUE){
        edges$color = factor(ifelse(edges$weight > 0, 'Positive', 'Negative'), levels=c('Positive', 'Negative'))
    }
    
    # select nodes by degree
    # ----------------------
    e = edges[edges$weight > 0,]
    g = graph.data.frame(e)
    if(!is.na(min_degree)){
	v.k = names(which(degree(g) >= min_degree))
	v.n = unlist(sapply(v.k, function(a) names(neighbors(g, a))))
	v.k = sort(unique(c(v.k, v.n)))
	v.r = setdiff(names(V(g)), v.k)
	g = delete_vertices(g, v.r)
    }
    n = V(g)$name
    
    # calculate graph layout (positive nodes)
    # ---------------------------------------
    e = edges[edges$weight > 0,]
    e = e[(e$source %in% n) & (e$target %in% n),]
    g = graph.data.frame(e)
    if(layout == 'fruchtermanreingold'){
        l = layout_with_fr(g, weights=e$weight)
    } else if(layout == 'circle'){
        l = layout_in_circle(g)
    } else if(layout == 'nicely'){
        l = layout_nicely(g)
    } else if(layout == 'bipartite'){
        l = layout_as_bipartite(g)
    } else if(layout == 'gem'){
        l = layout_with_gem(g)
    } else if(layout == 'mds'){
        l = layout_with_mds(g)
    } else {
        quit('error: invalid layout')
    }
    
    n = V(g)$name
    rownames(l) = n
    
    # calculate final graph (all nodes)
    # ---------------------------------
    e = edges[(edges$source %in% n) & (edges$target %in% n),]
    if(color_sign){e$color = factor(ifelse(e$weight > 0, 'Positive', 'Negative'), levels=c('Positive', 'Negative'))}
    e$weight = abs(e$weight)
    g = graph.data.frame(e)
    node_order = V(g)$name
    V(g)$node_size = node_sizes[node_order]
    V(g)$node_color = as.character(node_colors[node_order])
    V(g)$node_stroke = as.numeric(node_strokes[node_order])
    l = l[node_order,]
    d = ggnetwork(g, arrow.gap=arrow.gap, layout=l)
    d$vertex.names = d$name
    
    # scale sizes (after layout)
    # --------------------------
    i = (d$x == d$xend & d$y == d$yend)
    #d[!i, 'weight'] = rescale_vector(d[!i, 'weight'], edge_scale)
    if(is.null(node_limits)){
        node_limits = c(min(V(g)$node_size), max(V(g)$node_size))
    }
    
    # plot parameters
    # ---------------
    if(!is.null(lab.use)){
        d$vertex.names = ifelse(d$vertex.names %in% lab.use, as.character(d$vertex.names), '')
	d$label = ifelse(d$label %in% lab.use, as.character(d$label), '')
    }
    d$node_color = factor(d$node_color, levels=levels(node_colors))
    d$color = factor(d$color, levels=levels(edges$color))
    d$linetype = factor(d$linetype, levels=levels(edges$linetype))
    if(is.null(node_pal)){
        if(nlevels(d$node_color) == 1){node_pal = 'black'} else {node_pal = set.colors}
    }
    if(is.null(edge_pal)){
        if(nlevels(d$color) == 1){edge_pal = 'grey'} else {edge_pal = set.colors}
    }

    # self edges
    # ----------
    if(self == TRUE){
        self_nodes = edges[edges$source == edges$target, 'source']
	d$node_stroke = ifelse(i & d$source %in% self_nodes, d$node_stroke, 0)
    }

    if(color_sign){
        d = d[order(d$color, decreasing=TRUE),]
	d$alpha = ifelse(d$color == 'Positive', edge_alpha, .15*edge_alpha)
    } else {
        d$alpha = edge_alpha
    }
    
    # plot network
    # ------------
    p = ggplot(d)
    
    # scale nodes or edges (can't do both)
    if(size_legend == 'nodes'){
        p = p + geom_segment(data=d[!i,], aes(x=x, y=y, xend=xend, yend=yend, colour=color, linetype=linetype), size=d[!i,'weight'], alpha=d[!i,'alpha'], arrow=arrow) +
	    geom_point(data=d[i,], pch=21, aes(x=x, y=y, fill=node_color, size=node_size, stroke=node_stroke), alpha=node_alpha) +
	    scale_radius(size_title, limits=node_limits, range=node_scale, guide=size_guide, trans=size_trans)
    } else {
        p = p + geom_segment(data=d[!i,], aes(x=x, y=y, xend=xend, yend=yend, colour=color, linetype=linetype, size=weight), alpha=d[!i,'alpha'], arrow=arrow) +
	    geom_point(data=d[i,], pch=21, aes(x=x, y=y, fill=node_color, stroke=node_stroke), size=max(node_scale), alpha=node_alpha) +
	    scale_size(size_title, limits=edge_limits, range=edge_scale, guide=size_guide, trans=size_trans)
    }

    # construct the rest of the plot

    # **new** fix overlapping node/edge labels
    lab = d; lab$pad = NA; lab$size = NA
    lab[i,]$label = lab[i,]$vertex.names
    lab[i,]$pad = node_padding
    lab[i,]$size = node_lab.size
    lab[!i,]$x = (lab[!i,]$x + lab[!i,]$xend)/2
    lab[!i,]$y = (lab[!i,]$y + lab[!i,]$yend)/2
    lab[!i,]$pad = edge_padding
    lab[!i,]$size = edge_lab.size
        
    p = p + 
	#geom_text_repel(data=d[i,], aes(x=x, y=y, label=vertex.names), box.padding=unit(node_padding, 'lines'), segment.color=segment.color, size=node_lab.size) +
	#geom_text_repel(data=d[!i,], aes(x=(x+xend)/2, y=(y+yend)/2, label=label), box.padding=unit(edge_padding, 'lines'), segment.color=segment.color, size=edge_lab.size) +
	geom_text_repel(data=lab, aes(x=x, y=y, label=label), box.padding=unit(max(node_padding, edge_padding), 'lines'), segment.color=segment.color, size=lab$size, max.overlaps=max.overlaps) +
	scale_colour_manual(edge_title, breaks=levels(d$color), labels=levels(d$color), values=edge_pal, guide=edge_guide, drop=drop_levels) +
	scale_linetype_manual(linetype_title, guide=linetype_guide, drop=drop_levels, values=linewidth_pal) +
	scale_fill_manual(node_title, breaks=levels(d$node_color), labels=levels(d$node_color), values=node_pal, guide=node_guide, drop=drop_levels) +	
	theme_blank() + ggtitle(ggtitle)
    
    # fix legends
    # -----------
    p = p + guides(
                fill = guide_legend(override.aes = list(size=.75*max(node_scale))),
		colour = guide_legend(override.aes = list(size=2))
		)
    if(do.legend == FALSE){
        p = p + theme(legend.position='none')
    }
    if(nlevels(d[i, 'node_color']) == 1){
        p = p + guides(fill=FALSE)
    }
    if(nlevels(d[!i, 'color']) == 1){
        p = p + guides(colour=FALSE)
    }
    if(size_legend == 'nodes'){
        if(min(d[i, 'node_size'], na.rm=T) == max(d[i, 'node_size'], na.rm=T)){p = p + guides(size=FALSE)}
    }
    if(nlevels(d[!i, 'linetype']) == 1){
        p = p + guides(linetype=FALSE)
    }
    p
    
    # write output
    # ------------
    if(!is.null(out)){
        save_plot(p, file=out, nrow=nrow, ncol=ncol)
    }
    if(ret.data == FALSE){
        p        
    } else {
        list(p=p, d=d)
    }
}
