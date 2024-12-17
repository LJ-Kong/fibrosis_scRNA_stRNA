

load_maps = function(h='human_genes'){

    # Get gene lists
    hgenes = readLines(get_path('code/human.genes.txt'))
    
    # Get gene synonyms
    hsyn = get_path('code/human.gene_map.txt')
    hsyn = read.table(hsyn, row.names=1, sep='\t', stringsAsFactors=F, quote='', comment.char='')

    # Get other symbols
    hsym = get_path('code/human.db2sym.txt')
    hsym = read.table(hsym, row.names=1, sep='\t', stringsAsFactors=F, quote='', comment.char='')
    
    # Combine synonyms and symbols
    hsyn = unique(rbind(hsyn, hsym))
    
    # Return maps
    return(list(h=h, hgenes=hgenes, hsyn=hsyn, hsym=hsym))
}


maps = load_maps()
list2env(maps, .GlobalEnv)

fix_names = function(names){
    names = toupper(names)
    names = gsub('[^a-zA-Z0-9]', '', names)
    return(names)
}

get_synonyms = function(genes, target='human', do.unlist=TRUE){
    genes = fix_names(genes)
    if(target == 'human'){genes = hsyn[genes,1]}
    if(do.unlist == TRUE){
        genes = unlist(strsplit(genes, ','))
    } else {
        genes = strsplit(genes, ',')
    }
    genes
}

map_gene = function(genes, target='human', source='human', do.unlist=TRUE){
    get_synonyms(genes, target=target, do.unlist=do.unlist)
}


get_functions = function(genes, type='desc', org='human'){
    res = setNames(rep('', length(genes)), genes)
    if(type == 'desc'){
        if(org == 'human'){
            hmart = useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
	    anno = getBM(attributes=c('hgnc_symbol', 'description'), mart=hmart)
    	}
    	anno = anno[anno[,1] != '',]
    	anno = anno[!duplicated(anno[,1]),]
    	anno = data.frame(anno, row.names=1)
    	anno[,1] = gsub('\\[.*\\]$', '', anno[,1])
	i = intersect(rownames(anno), genes)
	res[i] = anno[i,1]
    }
    if(type == 'summary'){
        if(org == 'human'){
	    entrez = mapIds(org.Hs.eg.db, genes, 'ENTREZID', 'SYMBOL')
	    res[] = unlist(getGenes(entrez, fields='summary')$summary)
	}
    }
    return(res)
}
