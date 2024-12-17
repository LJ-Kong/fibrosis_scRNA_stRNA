source('../load_data.r')

# ------------------------
# 1. write metadata tables
# ------------------------

cmeta = as.data.table(vis$meta.data)
cmeta = cmeta[, .(Barcode=gsub('.*\\.', '', rownames(vis$meta.data)), "Shorthand ID"=name, Subject_PubID=pubid, Diagnosis=dx, nGene=nGene, nUMI=nUMI, Region=region, Cluster=cluster)]
write.table(cmeta, file='../table/Table_S6.spot_meta.tsv', sep='\t', quote=F, row.names=F)


smeta = as.data.table(vis$meta.data)
smeta = smeta[,num_spot := .N,name][,lapply(.SD, meanmode),name][,.("Visium Slide ID"=sample, "Disease Label"=dx, "Patient ID"=pubid, "Shorthand ID"=name, Inflamed=inflamed, Spots=num_spot, "Mean nGene"=nGene, "Mean nUMI"=nUMI)] 

i = vis$meta.data$dx == 'Adjacent'
all_stats = c('Total', '', '', '', ncol(vis$counts), mean(colSums(vis$counts > 0)), mean(colSums(vis$counts)))
non_stats = c('Non-stricturing', '', '', '', ncol(vis$counts[,i]), mean(colSums(vis$counts[,i] > 0)), mean(colSums(vis$counts[,i])))
fib_stats = c('Stricturing', '', '', '', ncol(vis$counts[,!i]), mean(colSums(vis$counts[,!i] > 0)), mean(colSums(vis$counts[,!i])))
smeta = rbind(as.data.frame(smeta), rep('', ncol(smeta)), all_stats, non_stats, fib_stats)

write.table(smeta, file='../table/Table_S6.slide_meta.tsv', sep='\t', quote=F, row.names=F)


# ---------------------------------------
# 2. load and process additional datasets
# ---------------------------------------

# IBD GWAS gene sets
# ------------------
gwas.all = load_signature(get_path('code/gwas.txt'))
gwas.ibd = gwas.all$IBD
gwas.cd = gwas.all$CD
gwas.uc = gwas.all$UC
gwas.mono = gwas.all$mono

# single-cell metadata
# --------------------

# extract single-cell data
sco.ident = i2f[sco$meta.data$annotation2]
sco.sample = sco$meta.data$sample
sco.health = sco$meta.data$health
sco.status = sco$meta.data$status
sco.fraction = sco$meta.data$fraction
sco.procedure = sco$meta.data$procedure

# calculate biopsy/resection
epi = sco.fraction == 'epi'
imu = sco.fraction == 'imu'


# receptor-ligand interactions
# ----------------------------

# load ligand-receptor list (fantom5)
ligrec = fread(get_path('code/ligand_receptor.literature.txt'))
colnames(ligrec) = c('lig', 'rec')

# subset to highly expressed genes
gep.norm = gep.norm[, apply(gep.norm, 2, max) >= 1]
ligrec = ligrec[apply(ligrec, 1, function(a) all(a %in% colnames(gep.norm))),]
lig = ligrec$lig
rec = ligrec$rec

# permute ligand-receptor pairs
ligrec.shuf = sapply(1:100, function(a){a = ligrec; a$lig = sample(a$lig, replace=F); a$rec = sample(a$rec, replace=F); a}, simplify=F)

stop()

# ------------------------------------
# 3. cluster and annotate visium spots
# ------------------------------------


# cluster spots by gene expression
# --------------------------------

if(FALSE){# set to TRUE to run visium spot clustering pipeline

# run visium spot clustering on counts matrix
vis = run_seurat(counts=vis$counts, name='test', num_genes=2000, num_pcs=25, do.batch='liger', write_out=F, verbose=T)

# cluster visium spots based on PCs 1-25
v = run_graph_cluster(data=vis$pca.rot[,1:25], k=250, do.fast=T, method='louvain')
vis$meta.data$liger.nn1.louvain.cosine.k250 = v[,1]

}


# euclidean distances between cell types
# --------------------------------------

get_proximity = function(v, min_spots, type='mean'){

    # map spots to samples
    spot2sample = table(v[,1], vis$meta.data$sample)
    
    # calculate distances
    d = vis$impos
    V = naturalsort(unique(v[,1]))
    D = sapply(names(d), function(a){
        b = d[[a]]
        b = b[intersect(rownames(b), rownames(vis$meta.data)),]
        g = v[rownames(b),]
        if(sum(table(g)) > 0){
            b = as.matrix(dist(b))
            r = nice_agg(t(nice_agg(b, g, type)), g, type)
        } else {
            r = matrix(NA, nrow=length(V), ncol=length(V))
        }
        R = matrix(NA, nrow=length(V), ncol=length(V))
        rownames(R) = colnames(R) = V
        R[rownames(r),colnames(r)] = as.matrix(r)
        spots.rmv = names(which(spot2sample[,a] < min_spots))
        R[spots.rmv, spots.rmv] = NA
        R
    }, simplify=F)
    D.u = sapply(1:ncol(D[[1]]), function(j) rowMeans(sapply(D, function(di) di[,j]), na.rm=T))
    D.v = sapply(1:ncol(D[[1]]), function(j) apply(sapply(D, function(di) di[,j]), 1, var, na.rm=T))
    h = hclust(as.dist(D.u))
    return(list(D=D, D.u=D.u, D.v=D.v, h=h))
}


# annotate tissue regions based on cell-cell proximities
# ------------------------------------------------------

# get cell-cell proximities
v = vis$meta.data[,'liger.nn1.louvain.cosine.k250',drop=F]
d = get_proximity(v, min_spots=25)

# annotate tissue regions (epithelial/lp/immune/muscle)
spot2loc = factor(cutree(d$h, 4))
levels(spot2loc) = c('LP', 'Muscle', 'Epithelial', 'Immune')
spot2loc = factor(spot2loc, levels=c('Epithelial', 'LP', 'Muscle', 'Immune'))
vis$meta.data$region = spot2loc[vis$meta.data$liger.nn1.louvain.cosine.k250]


# annotate spot clusters based on differentially expressed genes
# --------------------------------------------------------------

# calculate differentially expressed genes
markers = p.find_all_markers(vis, ident.use=vis$meta.data$liger.nn1.louvain.cosine.k250, test.use='mast', min_fc=2)
saveRDS(markers, file='visium_markers.rds')

# annotate spots based on cells + gene expression
imap = list(
    '4'='Mature',
    '5'='TA',
    '6'='Stem',
    '1'='Adipose',
    '2'='Myeloid',
    '9'='Glia',
    '10'='Blood vessel',
    '11'='Neuron',
    '12'='Endothelium',
    '15'='Plasma cell',
    '3'='Inflammatory',
    '7'='ICC',
    '8'='Myenteric plexus',
    '14'='Muscle',
    '13'='Follicle'
)

# save spot annotations
vis$meta.data$cluster = factor(vis$meta.data$liger.nn1.louvain.cosine.k250, levels=names(imap))
levels(vis$meta.data$cluster) = unlist(imap[levels(vis$meta.data$cluster)])

# colorblind friendly plots (from dittoSeq)
# -----------------------------------------
CB=T


# ---------------------------------------------------------------------------------
# FIGURE 4: plot samples: h&e, deconvolution, spot clusters, and changes in disease
# ---------------------------------------------------------------------------------


# Figure 4A: plot representative samples (h&e, clusters, regions, deconvolution)
# ------------------------------------------------------------------------------

# initialize variables
sh = 'V10A14-143_A'
sd = 'V10A14-143_B'

se = 'V10A14-143_A'
si = 'V10A14-143_D'
sf = 'V11Y24-011_D'

cells.use = c('Epithelial', 'Fibroblast', 'Myeloid', 'B cell', 'T cell', 'Pericyte', 'Glia')

# plot deconvolution and clusters for select samples
scores = as.data.frame(q[,cells.use])
scores$Cluster = vis$meta.data$cluster
scores$Region = vis$meta.data$region
scores = scores[,c('Cluster', 'Region', setdiff(colnames(scores), c('Cluster', 'Region')))]
p1 = plot_tsne(vis, ident=F, coords=vis$impos[[sh]], scores=scores, do.sort=T, num_col=10, theme_void=T, do.label=F, add_image=vis$image[[sh]], image_ds=1, cblind=CB)
p2 = plot_tsne(vis, ident=F, coords=vis$impos[[sd]], scores=scores, do.sort=T, num_col=10, theme_void=T, do.label=F, add_image=vis$image[[sd]], image_ds=1, cblind=CB)
ps = plot_grid(p1,p2,nrow=2)
save_plot(ps, file='Fig_4A.repr_deconvolution.pdf', nrow=1.5, ncol=4)


# Figure 4B: plot all samples (clusters and regions)
# --------------------------------------------------

p1 = plot_images(vis, scores=data.frame(Cluster=vis$meta.data$cluster), do.image=F, ident=F, ncol=10, title.size=10, add_image=FALSE, cblind=CB)
save_plot(p1, file='Fig_4B.spot_clusters.pdf', nrow=1.5, ncol=4)

p2 = plot_images(vis, do.image=T, ident=F, ncol=10, title.size=10, cblind=CB)
save_plot(p2, file='Fig_4B.h&e_images.pdf', nrow=1.5, ncol=4)

# ------------------------------------------
# Table S7: validate cell type deconvolution
# ------------------------------------------
mi = sapply(colnames(p), function(a) p.find_markers(vis, a, ident.use=factor(ifelse(1:ncol(vis$data) %in% order(-p[,a])[1:1000], a, 'Other'), levels=c('Other', a)), test.use='mast', min_fc=2.5, dir='pos'), simplify=F)
saveRDS(mi, file='deconv_markers.rds')
write.table(t(sapply(mi, function(a) a[padjD < .05][alpha > .05][order(-log2fc)][1:25]$gene)), file='../table/Table_S7.deconv_markers.sort_fc.tsv', sep='\t', quote=F, col.names=F)
write.table(t(sapply(mi, function(a) a[padjD < .05][alpha > .05][order(padjD)][1:25]$gene)), file='../table/Table_S7.deconv_markers.sort_p.tsv', sep='\t', quote=F, col.names=F)

# Figure 4C: cell type proportions of spots
# -----------------------------------------

cluster2region = data.frame(unique(data.frame(vis$meta.data[!is.na(vis$meta.data$cluster),c('cluster', 'region')])), row.names=1)
w = nice_agg(p, vis$meta.data$cluster, 'mean')
p1 = simple_scatter(rep(1, nrow(w)), rev(1:nrow(w)), col=factor(levels(vis$meta.data$cluster), levels=levels(vis$meta.data$cluster)), pal=set.colors) + theme_void() + theme(legend.position='none')
p2 = ggheatmap(scale(w[,apply(w, 2, max) >= .01]), Rowv='none')

# Figure 4D: select marker genes for spot clusters
# ------------------------------------------------

genes.use = c('ADIPOQ', 'CIDEA', 'LPL', 'RBP4', 'PLIN1',
              'HBA2', 'HBB', 'SLC7A2', 'RERGL', 'GJA5',
	      'UCHL1', 'ELAVL4', 'NMU', 'SST', 'HTR3A',
	      'ANO1', 'KIT', 'BCAN', 'NKX3-2', 'LY6H',
	      'GRP', 'SMPX', 'HAND1', 'HAND2', 'TUBB2A',
	      'SPIB', 'LTB', 'TCL1A', 'TIGIT', 'TNFRSF13C'
	      )
scores = score_cells(vis, genes.use, groups=vis$meta.data$cluster)
p3 = ggheatmap(scale(scores), Rowv='none', Colv='none') + theme(axis.text.y=element_blank())

# Figure 4E: changes in cluster proportions with disease
# ------------------------------------------------------

# Calculate spot proportions per sample [samples x spot clusters]
u = t(table(vis$meta.data$cluster, vis$meta.data$sample))
u = u/rowSums(u)

# Get disease annotations for each sample
v = data.frame(unique(vis$meta.data[,c('sample', 'dx')]), row.names=1)
v = ifelse(v[rownames(u),,drop=F][,1] == 'Adjacent', 'Adjacent', 'Stricture')

# Aggregate spot proportions by disease label
d = as.data.frame(t(nice_agg(as.data.frame.matrix(u), v, 'mean'))) %>% rownames_to_column('Cluster')
d$Cluster = factor(d$Cluster, levels=rev(d$Cluster))
D = d %>% gather(Type, Frequency, -Cluster)

# Plot barplot showing spot proportions in health and disease
p4 = ggplot(d) + geom_segment(aes(x=0, xend=pmin(Adjacent, Stricture), y=Cluster, yend=Cluster), color='grey', size=1) +
            geom_segment(aes(x=pmin(Adjacent, Stricture), xend=pmax(Adjacent, Stricture), y=Cluster, yend=Cluster, color=ifelse(Adjacent > Stricture, 'Adjacent', 'Stricture')), size=1) + 
	    geom_point(data=D, aes(x=Frequency, y=Cluster, color=Type), size=1.5) + 
	    theme_cowplot() + xlab('') + ylab('') + scale_color_manual('', values=set.colors[c(2,1)]) + theme(axis.text.y=element_blank()) + theme(legend.position='none')

# Save Figures 4C, 4D, and 4E
# ---------------------------
ps = plot_grid(p1,p2,p3,p4,nrow=1,rel_widths=c(.2,1,1,.5),align='h')
save_plot(ps, file='Fig_4CDE.spot_clusters_composition.pdf', nrow=1, ncol=1.75)


# --------------------------------------------------
# 5. Extended Data Figure 5: collagen-hi fibroblasts
# --------------------------------------------------


# Figure S5A: comparison of inflammatory and collagen-hi fibroblasts
# ------------------------------------------------------------------

# differential expression tests to identify unique marker genes
sco.group = f2g[sco.ident]
group.inf = group.col = ifelse(sco.group == 'Fibroblast', 'Fibroblast', sco.ident)
group.inf[sco.ident == 'Inflammatory fibroblasts'] = 'Inflammatory fibroblasts'
group.col[sco.ident == 'Collagen-hi tissue fibroblast'] = 'Collagen-hi tissue fibroblast'
m.inf = p.find_markers(sco, 'Inflammatory fibroblasts', 'Fibroblast', test.use='mast', dir='both', ident.use=group.inf)
m.col = p.find_markers(sco, 'Collagen-hi tissue fibroblast', 'Fibroblast', test.use='mast', dir='both', ident.use=group.col)

# volcano plots of inflammatory and collagen-hi fibroblast genes
p1 = simple_scatter(m.inf$coefD, -log(m.inf$padjD), lab=m.inf$gene, lab.n=25, xlab='Fold change\n(Inflammatory vs. Other)', ylab='-log(adjusted P)', lab.size=3, xskip=c(-1,1)) + theme(legend.position='none')
p2 = simple_scatter(m.col$coefD, -log(m.col$padjD), lab=m.col$gene, lab.n=25, xlab='Fold change\n(Inflammatory vs. Other)', ylab='-log(adjusted P)', lab.size=3, xskip=c(-1,1)) + theme(legend.position='none')
ps = plot_grid(p1,p2,nrow=1)
save_plot(ps, file='Fig_S5AB.fibroblast_volcano.pdf', nrow=1.5, ncol=2)


# Figure S5B: median number of UMIs per cell type
# -----------------------------------------------

# calculate median number of UMIs per cell type
u = sort(tapply(sco$meta.data$nUMI, sco.ident, median), dec=T)

# compact visualization (display top 20 cell types)
i = f2g[names(u)] == 'Fibroblast'
a = names(u[!i][1:(19-sum(i))])
b = names(which(i))
u = names(u[names(u) %in% c(a,b)])
d = data.frame(ident=sco.ident, nUMI=sco$meta.data$nUMI)
d[! d$ident %in% u, 'ident'] = 'Other'
d$ident = factor(d$ident, levels=c(u, 'Other'))

# barplot showing median number of UMIs per cell type
p1 = ggplot(d) + geom_boxplot(aes(x=ident, y=log10(nUMI), fill=f2g[as.character(ident)])) + theme_cowplot() + theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) + xlab('') + scale_fill_manual('', values=set.colors)
save_plot(p1, file='Fig_S5C.numi_boxplot.pdf', nrow=1.5, ncol=2)



# -------------------------------------------------------------
# FIGURE 5A-C: Radial axis and spatial clustering of cell types
# -------------------------------------------------------------


# Estimate radial axis from epithelium to muscle layer
# ----------------------------------------------------

# first calculate line segments from epithelial to muscle spots
d = do.call(rbind,
    sapply(names(vis$impos), function(a){
        i = grep(a, rownames(vis$tsne.rot))
	b = nice_agg2(vis$tsne.rot[i,], vis$meta.data[i,]$region, FUN=function(a){geom_mean(replace_na(a, mean(a)))})
	d = data.frame(x0=b[1,1], x1=b[3,1], y0=b[1,2], y1=b[3,2])
	d$m=(d$y1-d$y0)/(d$x1-d$x0)
	d$b=d$y0-d$m*d$x0
	d
    }, simplify=F)
)

# next, project all spots along these radial axes
# proj$d1 = spatial axis projection
proj = sapply(unique(vis$meta.data$sample), function(a){
    i = vis$meta.data$sample == a
    b = scale(t(project_points(vis$tsne.rot[i,], slope=d[a,'m'], intercept=d[a,'b'])))
    u = tapply(b[,'d1'], vis$meta.data[i, 'region'], mean)
    if(u[['Epithelial']] > u[['Muscle']]){b[,'d1'] = -1*b[,'d1']}
    b
}, simplify=F)
proj = as.data.frame(do.call(rbind, proj)[colnames(vis$data),])

# compare to other "radial axis" definitions
# ------------------------------------------
epi = rownames(vis$meta.data)[vis$meta.data$cluster %in% c('Stem', 'TA', 'Mature')]
mus = rownames(vis$meta.data)[vis$meta.data$cluster %in% c('Muscle')]
di = sapply(vis$impos, function(a) as.matrix(dist(a)), simplify=F)
get_radial_axis = function(q){
    epi.dist = Reduce(c, sapply(di, function(a) setNames(apply(a[, colnames(a) %in% epi], 1, quantile, q), rownames(a)), simplify=F))[rownames(vis$meta.data)]
    mus.dist = Reduce(c, sapply(di, function(a) setNames(apply(a[, colnames(a) %in% mus], 1, quantile, q), rownames(a)), simplify=F))[rownames(vis$meta.data)]
    epi.dist/(epi.dist + mus.dist)
}
axis.0 = get_radial_axis(0)[colnames(vis$data)]
axis.5 = get_radial_axis(0.5)[colnames(vis$data)]

i = vis$meta.data$dx == 'Adjacent'
sh = sort(unique(vis$meta.data[i, 'sample']))
sd = sort(unique(vis$meta.data[!i, 'sample']))

scores = data.frame("radial_axis"=proj$d1, "new_definition_min"=axis.0, "new_definition_median"=axis.5, "EPCAM"=vis$data['EPCAM',], "ACTA2"=vis$data['ACTA2',])
p1 = plot_images(vis, scores=scores, num_col=6, image.use=c(sh[c(1,2,4)], sd[c(2,3,4)]), do.sort=T, pt.size=.5, ident=F, ncol=1, cblind=CB)
save_plot(p1, file='RFig.radial_axis.png', nrow=2.5, ncol=2.5)


# Figure 5A: spatial distribution of cell types
# ---------------------------------------------

# aggregate epithelial and muscle cells for plots
epi = setdiff(names(which(f2g == 'Epithelial')), c('Epithelial', 'Epithelium'))
mus = grep('[Pp]ericyte|[Gg]lia|SIP|Myofib', colnames(p), value=T)
P = as.data.frame(p[,setdiff(colnames(p), c(epi, mus))])
P$Epithelium = rowSums(p[,epi])
P$Muscle = rowSums(p[,mus])

# calculate radial axis position and moran's I for each cell type
d = as_adj(nng(vis$tsne.rot, k=25, use.fnn=T))
u = sapply(as.data.frame(P), function(a) mean(a*proj$d1))
v = sapply(as.data.frame(P), function(a) fast_moran_i(a, d))
w = p.adjust(apply(P, 2, function(a) tryCatch({wilcox.test(proj$d1[a > .10], proj$d1[sample(a, replace=F) > .10])$p.value}, error=function(e){1})), 'fdr')
f2g[['Epithelium']] = 'Epithelial'
f2g[['Muscle']] = 'Muscle'
i = names(which(w < .05))

p1 = simple_scatter(u[i], v[i], lab=i, lab.n=25, col=f2g[i], pal=set.colors, size=2, groups=1, xlab='Epithelium-Muscle axis (scaled)', ylab="Moran's I", xmin.n=2, xmax.n=2)
save_plot(p1, file='Fig_5A.spatial_stats.pdf', nrow=2, ncol=2)

# compare to alternative radial axis
U = sapply(as.data.frame(P), function(a) mean(a*scale(axis.5)))
p1 = simple_scatter(u[i], U[i], xmin.n=2, xmax.n=2, ymin.n=2, ymax.n=2, lab=i, lab.n=10, lab.type=c('up', 'down'), lm=T, xlab='Original definition', ylab='Alternative definition (Median)') + theme(legend.position='none')
save_plot(p1, file='Fig_R.radial_axis.scatter.png', nrow=2.5, ncol=2.5)


# Figure S4: spatial distributions of cell types
# ----------------------------------------------

i = vis$meta.data$dx == 'Adjacent'
sh = sort(unique(vis$meta.data[i, 'sample']))
sd = sort(unique(vis$meta.data[!i, 'sample']))
pt_size = .1

# plot clusters, radial axis, deconvolution for all tissue sections (non-stricturing)
cols.use = c('B cell', 'Epithelial', 'Fibroblast', 'Myeloid', 'Pericyte', 'T cell')
p1 = plot_images(vis, do.image=T, ident=F, ncol=10, do.title=F, image.use=sh, do.legend=F, cblind=CB)
p2 = plot_images(vis, scores=data.frame(Cluster=vis$meta.data$cluster), do.image=F, ident=F, ncol=10, do.title=F, add_image=F, image.use=sh, pt.size=pt_size, do.legend=F, cblind=CB)
p3 = plot_images(vis, scores=data.frame("Radial axis"=proj$d1), do.image=F, ident=F, ncol=10, do.title=F, add_image=F, image.use=sh, pt.size=pt_size, vmin=min(proj$d1), vmax=max(proj$d1), do.legend=F, cblind=CB)
ps = sapply(cols.use, function(a) {
        vmin = quantile(q[,a], .01)
	vmax = quantile(q[,a], .99)
	print(c(a, vmin, vmax))
        plot_images(vis, scores=q[,a], do.image=F, ident=F, ncol=10, do.title=F, do.sort=T, add_image=F, image.use=sh, pt.size=pt_size, vmin=vmin, vmax=vmax, do.legend=F, cblind=CB)
     }, simplify=F)
ps1 = plot_grid(p1, p2, p3, ps[[1]], ps[[2]], ps[[3]], ps[[4]], ps[[5]], ps[[6]], align='hv', ncol=1)

# plot clusters, radial axis, deconvolution for all tissue sections (non-stricturing)
cols.use = c('B cell', 'Epithelial', 'Fibroblast', 'Myeloid', 'Pericyte', 'T cell')
p1 = plot_images(vis, do.image=T, ident=F, ncol=10, do.title=F, image.use=sd, do.legend=F, cblind=CB)
p2 = plot_images(vis, scores=data.frame(Cluster=vis$meta.data$cluster), do.image=F, ident=F, ncol=10, do.title=F, add_image=F, image.use=sd, pt.size=pt_size, do.legend=F, cblind=CB)
p3 = plot_images(vis, scores=data.frame("Radial axis"=proj$d1), do.image=F, ident=F, ncol=10, do.title=F, add_image=F, image.use=sd, pt.size=pt_size, vmin=min(proj$d1), vmax=max(proj$d1), do.legend=F, cblind=CB)
ps = sapply(cols.use, function(a) {
        vmin = quantile(q[,a], .01)
	vmax = quantile(q[,a], .99)
	print(c(a, vmin, vmax))
        plot_images(vis, scores=q[,a], do.image=F, ident=F, ncol=10, do.title=F, do.sort=T, add_image=F, image.use=sd, pt.size=pt_size, vmin=vmin, vmax=vmax, do.legend=F, cblind=CB)
     }, simplify=F)
ps2 = plot_grid(p1, p2, p3, ps[[1]], ps[[2]], ps[[3]], ps[[4]], ps[[5]], ps[[6]], align='hv', ncol=1)

ps = plot_grid(ps1, ps2, nrow=2)
ggsave('Fig_S4.all_tissue_sections.pdf', plot=ps, width=6.5, height=9, bg='white', units='in')


# Figure 5B: examples of spatially clustered cell types
# -----------------------------------------------------

cells.use = c(
    'Inflammatory fibroblasts'='V10A14-143_D',
    'Collagen-hi tissue fibroblast'='V10S15-054_B',
    'CD63+CD81+ macrophage'='V10A14-143_C'
    )

# plot cell type density and radial axis
p1 = plot_tsne(vis, scores=data.frame(ident=p[,names(cells.use)[[1]]]), do.sort=T, theme_void=T, s.use=cells.use[[1]], ident=F, num_col=1, vmax=0.7, cblind=CB)
p2 = plot_tsne(vis, scores=data.frame(ident=p[,names(cells.use)[[2]]]), do.sort=T, theme_void=T, s.use=cells.use[[2]], ident=F, num_col=1, vmax=0.7, cblind=CB)
p3 = plot_tsne(vis, scores=data.frame(ident=p[,names(cells.use)[[3]]]), do.sort=T, theme_void=T, s.use=cells.use[[3]], ident=F, num_col=1, vmax=0.7, cblind=CB)
p4 = plot_tsne(vis, scores=data.frame(proj=proj$d1), do.sort=T, theme_void=T, s.use=cells.use[[1]], ident=F, num_col=1, vmin=-2, vmax=2, cblind=CB)
p5 = plot_tsne(vis, scores=data.frame(proj=proj$d1), do.sort=T, theme_void=T, s.use=cells.use[[2]], ident=F, num_col=1, vmin=-2, vmax=2, cblind=CB)
p6 = plot_tsne(vis, scores=data.frame(proj=proj$d1), do.sort=T, theme_void=T, s.use=cells.use[[3]], ident=F, num_col=1, vmin=-2, vmax=2, cblind=CB)
ps = plot_grid(plot_grid(p1, p4, nrow=2), plot_grid(p2, p5, nrow=2), plot_grid(p3, p6, nrow=2), nrow=1)
save_plot(ps, file='Fig_5B.radial_axis_examples.pdf', nrow=1.25, ncol=1.25)


# Figure 5C: distribution of radial axis positions
# ------------------------------------------------
cells.use = c(names(which(abs(scale(u[grep('Epi|Mus', names(u), invert=T)])[,1]) > 1.5)), 'Inflammatory fibroblasts', 'Epithelium', 'Muscle')

P = as.data.frame(P)
d = do.call(rbind, sapply(c(0.1), function(b) do.call(rbind, sapply(colnames(P), function(a) data.frame(ident=rep(a, sum(P[,a] > b)), proj=proj$d1[P[,a] > b], cutoff=rep(b, sum(P[,a] > b))), simplify=F)), simplify=F))
d$group = factor(ifelse(d$ident %in% cells.use, d$ident, 'Other'), levels=c(setdiff(names(sort(u[cells.use])), 'Other'), 'Other'))
d = d[d$group != 'Other',]
d$type = ifelse(d$group %in% names(which(u[cells.use] < 0)), 'Epithelial', 'Muscle')
p1 = ggplot(d) + geom_density(aes(x=proj, color=group), linewidth=1.25, adjust=1.5) + theme_cowplot() + scale_color_manual('Cell type', values=brewer.pal(11, 'Set3')) + xlab('Radial Axis') + ylab('Density') + facet_grid(type ~ .)
save_plot(p1, file='Fig_5C.radial_axis_density.pdf', nrow=1.25, ncol=1.25)

# get projection coordinates
cells.use = names(sort(tapply(d$proj, d$ident, mean)))
d = as.data.table(d)
d$ident = factor(d$ident, levels=cells.use)
pp = sapply(cells.use, function(a) wilcox.test(proj$d1[P[,a] > .1], proj$d1[P[,a] < .1])$p.value) # all tests significant
pp = data.frame(x=factor(names(pp), levels=cells.use), y=2.25, lab='***')
p2 = ggplot(d) + geom_boxplot(aes(x=ident, y=proj, fill=ident)) + xlab('Cell type') + ylab('Radial axis') + scale_fill_manual('', values=set.colors) + theme_cowplot() + geom_text(data=pp, aes(x=x, y=y, label=lab)) + theme(axis.text.x=element_blank())
save_plot(p2, file='Fig_S5D.radial_axis_boxplot.pdf', nrow=1.25, ncol=1.5)


# -----------------------------------------------------
# FIGURE 5D-G: Subclustering the enteric nervous system
# -----------------------------------------------------


# sub-cluster neuron spots and calculate marker genes
# ---------------------------------------------------

# identify neuron spots based on ens marker genes
m.neur = readLines(get_path('/code/ens_markers.txt'))
nscore = score_cells(vis, list(Neuron=m.neur), combine_genes='scale2')
neuron = rownames(nscore)[nscore[,1] > 1]    

if(FALSE){

    # sub-cluster neuron spots
    # ------------------------
    
    # load neurotransmitters and neuropeptides
    m.nts = readLines(get_path('code/neurosignaling_genes.txt'))
    g.use = intersect(rownames(vis$data), sort(unique(c(m.neur, m.nts))))
    
    # subcluster neurons using ens markers + neurosignaling genes
    neur = run_seurat(name='neur', counts=vis$counts[g.use, neuron], minc=0, ming=0, num_pcs=8, do.batch='combat', verbose=T, write_out=T)
    v.neur = run_graph_cluster(neur$pca.rot[,1:8], k=50, do.fast=T, method='infomap')
    ident.vis = ifelse(colnames(vis$data) %in% rownames(v.neur), v.neur[colnames(vis$data),1], NA)
    
    # calculate neuron markers
    #-------------------------
    
    m.neur = p.find_all_markers(neur, ident.use=v.neur[,1], test.use='mast', dir='both', min_fc=1.5)
    m.vis = p.find_all_markers(vis, ident.use=ident.vis, test.use='mast', dir='both', min_fc=1.5)
    mi = m.neur[padjH < 1e-5][alpha > .2][order(-log2fc)]
    mj = m.vis[padjH < 1e-5][alpha > .2][order(-log2fc)]
    Mi = tapply(mi$gene, mi$ident, function(a) a[1:min(length(a), 25)])
    Mj = tapply(mj$gene, mj$ident, function(a) a[1:min(length(a), 25)])
    
    # save neuron results
    res = list(neur=neur, v.neur=v.neur, m.neur=m.neur, m.vis=m.vis)
    saveRDS(res, file='neur.rds')

} else {

    # alternatively, load neuron sub-clustering results
    # -------------------------------------------------
    
    res = readRDS('neur.rds')
    neur = res$neur
    v.neur = res$v.neur
    m.neur = res$m.neur
    m.vis = res$m.vis
    ident.vis = ifelse(colnames(vis$data) %in% rownames(v.neur), v.neur[colnames(vis$data),1], NA)

}


# annotate neuron sub-clusters
# ----------------------------

# calculate neuron proximities
v = data.frame(ident.vis, row.names=rownames(vis$meta.data))
d = get_proximity(v, min_spots=5)

# annotate neuron sub-clusters
imap = c(
    '1'='ICC',
    '2'='Epithelial',
    '3'='Fibroblast',
    '4'='PMN', # PENK, CARTPT, CHGB, SLC5A7
    '5'='PSN 1', # SST, NMU, TAC1, HTR3A, GALR1
    '6'='PSVN', # SSTR1, SSTR2, NPY2R, VIP, NTSR1, GCGR
    '7'='Contamination',
    '8'='PSN 2', # NMU, NOG, SST, TAC1, HTR3A, GALR1
    '9'='EEC'
    )
i = ident.vis %in% names(imap)
ident.vis[i] = imap[ident.vis[i]]
ident.vis = factor(ident.vis, levels=c('PMN', 'PSN 1', 'PSN 2', 'PSVN', 'Epithelial', 'EEC', 'Fibroblast', 'ICC', 'Contamination'))
ident.vis = setNames(ident.vis, colnames(vis$data))

# Figure 5D - neuron spots colored by subcluster
# ----------------------------------------------
p1 = plot_tsne(neur, scores=data.frame("Neuron score"=nscore[colnames(neur$data),1], Cluster=ident.vis[colnames(neur$data)]), ymax=.99, num_col=3, do.label=F, cblind=CB)

# Figure 5E - representative plots of neuron subclusters
# ------------------------------------------------------
ident.new = factor(ifelse(grepl('^P', ident.vis), as.character(ident.vis), as.character(vis$meta.data$region)), levels=c('PMN', 'PSN 1', 'PSN 2', 'PSVN', 'Epithelial', 'LP', 'Muscle', 'Immune'))
p4 = plot_tsne(vis, scores=data.frame(Neuron=ident.new), ident=F, do.label=F, do.title=F, theme_void=T, dpal=c(set.colors[1:4], gray.colors(6)[3:6]), s.use='V10S15-055_A', na.value='#dddddd', do.revsort=T, pt.size=1.5, cblind=CB) + ggtitle('V10S15-055_A')

# Figure 5G - heatmap of marker genes for neuron subclusters
# ----------------------------------------------------------
p2 = ggdendrogram(d$h, rotate=90) + scale_y_reverse() + theme_void() + theme(axis.text.x=element_blank(), axis.text.y=element_blank())
genes.use = c('GRP', 'PENK', 'CARTPT', 'SLC5A7', 'NOS1', 'COLQ',
              'NMU', 'TAC1', 'SST', 'GALR1', 'HTR3A', 'CALCB', 'NOG',
	      'TH', 'VIP', 'GCGR', 'SSTR1', 'SSTR2', 'OPRM1',
	      'NTS', 'PYY', 'GCG', 'EPCAM', 'PIGR',
	      'COL1A1', 'COL1A2', 'RSPO3', 'GREM1',
	      'ANO1', 'KIT', 'GALR2', 'TACR2')
scores.use = nice_agg(as.matrix(t(2**vis$data[genes.use,]-1)), ident.vis, 'mean')
scores.use = log2(scores.use + 1)
Rowv = d$h$label[d$h$order]
p3 = ggheatmap(scale(scores.use), Rowv=Rowv)

# Save Figures 5D, 5E, and 5G
# ---------------------------
ps = plot_grid(p2,p3,nrow=1,rel_widths=c(.25, 1),align='hv')
ps = plot_grid(plot_grid(p1,p4,nrow=1,rel_widths=c(3,1)),plot_grid(p2,p3,rel_widths=c(.2,1)),nrow=2)
save_plot(ps, file='Fig5DEG.neuron_subclusters.pdf', nrow=3, ncol=3)

# Figure 5F - neuron-region enrichment (radial axis)
# --------------------------------------------------

# get projection coordinates
d = data.frame(ident=ident.vis, proj=proj$d1)
d = d[!is.na(d$ident),]
d = d[grep('^P', d$ident),]

# calculate p-values
pp = data.frame(t(sapply(c('PMN', 'PSN 1', 'PSN 2', 'PSVN'), function(a) c(x=a, y=2, pval=wilcox.test(d[d$ident == a,]$proj, d[d$ident != a,]$proj)$p.value))))
pp$pval = as.numeric(pp$pval)
pp$y = as.numeric(pp$y)
pp$lab = ifelse(pp$pval < .001, '***', ifelse(pp$pval < .01, '**', ifelse(pp$pval < .05, '*', 'N.S.')))

# plot spatial coordinates
p1 = ggplot(d) + geom_boxplot(aes(x=ident, y=proj, fill=ident)) + xlab('Neuron cluster') + ylab('Radial axis') + scale_fill_manual('', values=set.colors) + theme_cowplot() + theme(legend.position='none') + geom_text(data=pp, aes(x=x, y=y, label=lab))
save_plot(p1, file='Fig_5F.neuron_radial_axis.pdf', nrow=1.25, ncol=1)


# --------------------------------
# FIGURE 6: Cell-cell interactions
# --------------------------------


# load cell-cell interactions (sparcc)
# ------------------------------------

# load sparcc correlations
cci = readRDS('../2.SparCC/sparcc.res.rds')$Cor
cci[cci <= .4] = 0
rownames(cci) = colnames(cci) = colnames(p)

# select cell-cell interactions (highly correlated pairs)
cci.use = apply(which(cci > .4, arr.ind=T), 1, function(a) rownames(cci)[a])
j = apply(cci.use, 2, function(a) a[[1]] != a[[2]])
cci.use = cci.use[,j]


# Figure 6A - SparCC correlations
# -------------------------------

# calculate cell-cell frequencies
dx = vis$meta.data$dx
u = sapply(colnames(p), function(a) sapply(colnames(p), function(b) sum(p[dx == 'Adjacent',a]*p[dx == 'Adjacent',b])))
v = sapply(colnames(p), function(a) sapply(colnames(p), function(b) sum(p[dx != 'Adjacent',a]*p[dx != 'Adjacent',b])))
u = log(u/sum(u) + 1e-5)
v = log(v/sum(v) + 1e-5)

# plot cell-cell network
ncol = anno[,.(final, annotation3)]
ncol = setNames(ncol$annotation3, ncol$final)
ecol = ifelse(v - u > 1, 'Stricture', ifelse(v - u < -1, 'Adjacent', 'Other'))
eord = c('Adjacent', 'Stricture', 'Other')
epal = c(set.colors[c(2,1)], '#cccccc')

# select nodes to label
ecol[cci == 0] = 'Other'
lab1 = sort(unique(as.vector(apply(which(ecol == 'Stricture', arr.ind=T), 1, function(a) rownames(ecol)[a]))))
lab2 = c('Neutrophils', 'Enteric glia', 'Dendritic cell - cDC1', 'GZMB-hi B cell', 'Tissue fibroblast')
lab2 = names(which(apply(cci[lab2,], 2, any) | apply(cci[,lab2], 1, any)))
lab3 = c('Stem-like cells-OLFM4-hi', 'Paneth cells', 'Goblet cells-DEF-hi')
lab.use = c(lab1, lab2, lab3)

p1 = nice_network(graph=cci, node_colors=ncol, edge_colors=ecol, edge_colors_order=eord, edge_pal=epal, lab.use=lab.use)
save_plot(p1, file='Fig_6A.sparcc_network.pdf', nrow=2.5, ncol=2)


# Figure 6B - representative images of cell-cell interactions
# -----------------------------------------------------------

ix.use = list()
ix.use[[1]] = c('Stem-like cells-OLFM4-hi', 'Paneth cells', 'Goblet cells-DEF-hi')
ix.use[[2]] = c('Tissue fibroblast', 'Lyve1 macrophage', 'Venous endothelium')
ix.use[[3]] = c('Enteric glia', 'SIP syncytium fibroblast', 'Contractile pericytes')
ix.use[[4]] = c('B cells', 'Follicular helper T cells', 'CCR7hi CD4 T cells')
ix.use[[5]] = c('Monocyte-FPR1-hi', 'Inflammatory fibroblasts', 'Plasma cells-IgG')
ix.use[[6]] = c('CD63+CD81+ macrophage', 'CXCR6-hi CD4 T cells', 'Cycling T cells')

sa.use = list()

sa.use[[1]] = 'V10A14-143_A'
sa.use[[2]] = 'V10A14-143_A'
sa.use[[3]] = 'V10A14-143_A'
sa.use[[4]] = 'V10S15-054_D'
sa.use[[5]] = 'V10A14-143_D'
sa.use[[6]] = 'V10A14-143_C'

ps = sapply(1:length(ix.use), function(i){a=ix.use[[i]]
    plot_tsne(vis, multichannel=p[,a], theme_void=T, ident=F, coords=vis$impos[[sa.use[[i]]]], ymax=.995, pt.size=1, title.use=sa.use[[i]], cblind=CB)
}, simplify=F)

ps = plot_grid(plotlist=ps, align='hv')

save_plot(ps, file='Fig_6B.cci_examples.pdf', nrow=1.5, ncol=2.25)


# Figure 6C - statistics on cell-cell interactions
# ------------------------------------------------

# co-occurrence statistics
stats = as.data.table(data.frame(do.call(rbind, sapply(ix.use, function(b){t(sapply(unique(vis$meta.data$sample), function(a){i=vis$meta.data$sample == a; aa=mean(apply(p[i,b], 1, prod)); bb=prod(colMeans(p[i,b])); c(b, aa, bb, aa/bb, 'True')}))}, simplify=F))))
nulls = as.data.table(data.frame(do.call(rbind, sapply(1:10, function(ii) do.call(rbind, sapply(ix.use, function(b){t(sapply(unique(vis$meta.data$sample), function(a){i=vis$meta.data$sample == a; pb=p[i,b];
pb = apply(pb, 2, function(a) sample(a)); aa=mean(apply(pb, 1, prod)); bb=prod(colMeans(pb)); c(b, aa, bb, aa/bb, 'Null')}))}, simplify=F)), simplify=F))))
stats = rbind(stats, nulls)
colnames(stats) = c('a','b','c','p1','p2','fc','type')
stats$p1 = as.numeric(stats$p1)
stats$p2 = as.numeric(stats$p2)
stats$fc = as.numeric(stats$fc)
stats$type = factor(stats$type, levels=c('True', 'Null'))
stats[, ix := paste(a,b,c,sep='; ')]
names(ix.use) = sapply(ix.use, paste, collapse='; ')

# plot statistics
stats = stats[rev(1:nrow(stats)),]
stats[log(fc) > 4, 'fc'] = exp(4)
plab = data.frame(ix=names(ix.use), fc=exp(4.5), lab='***')
plab[grep('Infl', plab$ix), 'lab'] = '*'
p1 = ggplot(stats, aes(x=log(fc), y=factor(ix, levels=rev(names(ix.use))))) + geom_quasirandom(aes(color=type), groupOnX=FALSE) + geom_boxplot(data=stats[type == 'True'], width=.5, outlier.shape=NA, fill=NA) + theme_cowplot() + xlab('Colocalization logFC') + ylab('') + geom_vline(xintercept=0, linetype='dashed') +
geom_text(data=plab, aes(label=lab), angle=90) + scale_color_manual('', values=c('#000000', '#cccccc'))
save_plot(p1, file='Fig_6C.cci_colocalization.pdf', nrow=2, ncol=2)


# Calculate spatially clustered receptor-ligand interactions
# ----------------------------------------------------------
get_rli = function(c1, c2, niter){
    
    if(c1 == 'all'){
        p1 = rep(1, nrow(p))
	i1 = rep(T, length(lig))
	j1 = rep(T, length(rec))
    } else {
        p1 = p[,c1]
	i1 = (gep.norm[c1, lig] > 1) & (gep.norm[c1, rec] < 1)
	j1 = (gep.norm[c1, rec] > 1) & (gep.norm[c1, lig] < 1)
    }
    
    if(c2 == 'all'){
        p2 = rep(1, nrow(p))
	i2 = rep(T, length(rec))
	j2 = rep(T, length(lig))
    } else {
        p2 = p[,c2]
	i2 = (gep.norm[c2, rec] > 1) & (gep.norm[c2, lig] < 1)
	j2 = (gep.norm[c2, lig] > 1) & (gep.norm[c2, rec] < 1)
    }
    
    i = (i1 & i2)
    j = (j1 & j2)
    g1 = c(lig[i], rec[j])
    g2 = c(rec[i], lig[j])
    print(length(g1))

    if(length(g1) == 0){d = c(c1, c2, rep(NA, 10)); names(d) = c('c1', 'c2', 'true', 'shuf', 'sd', 'ratio', 'pval', 'epval', 'lig', 'rec', 'g1', 'g2'); return(t(d))}
    
    u = t(vis$data[g1,,drop=F])
    v = t(vis$data[g2,,drop=F])
    w = colSums(p1*u*p2*v)

    W = c()
    for(k in 1:niter){
        j = ave(1:ncol(vis$data), vis$meta.data$sample, FUN=sample)
	P2 = p2[j]
	V = v[j,,drop=F]
	W[[k]] = colSums(p1*u*P2*V)
    }
    W = do.call(cbind, W)
    d = data.table(data.frame(c1, c2, w, rowMeans(W), apply(W, 1, sd), w/rowMeans(W), rowSums(w < W)/ncol(W), dnorm(w, rowMeans(W), apply(W, 1, sd)), colnames(u), colnames(v), gep.norm[c1, g1], gep.norm[c2, g2]))
    colnames(d) = c('c1', 'c2', 'true', 'shuf', 'sd', 'ratio', 'pval', 'epval', 'lig', 'rec', 'g1', 'g2')
    d
}


# Figure 6D - identify significant receptor-ligand interactions
# -------------------------------------------------------------

# first calculate significant receptor-ligand interactions
if(FALSE){
    res = sapply(1:ncol(cci.use), function(a) get_rli(cci.use[1,a], cci.use[2,a], 100), simplify=F)
    res = do.call(rbind, res[sapply(res, typeof) == 'list'])
    res$epval = pnorm(res$true, mean=res$shuf, sd=res$sd, lower.tail=F)
    res[, epadj := p.adjust(epval, 'fdr')]
    res[, gene1 := lig]
    res[, gene2 := rec]
    res$lig = ifelse(res$gene1 %in% lig, res$gene1, res$gene2)
    res$rec = ifelse(res$gene1 %in% rec, res$gene1, res$gene2)
    res[, c1.alpha := ifelse(c1 < c2, c1, c2)]
    res[, c2.alpha := ifelse(c1 < c2, c2, c1)]
    f2g[['Neutrophils']] = 'Myeloid'
    res[, group1 := f2g[c1]]
    res[, group2 := f2g[c2]]
    saveRDS(res, file='cci_ligrec.rds')
} else {
    res = readRDS('cci_ligrec.rds')
}

# plot interaction network of follicle/inflammatory hits
cells.use = setdiff(names(which(apply(scale(nice_agg(p, vis$meta.data$cluster, 'mean'))[c('Follicle', 'Inflammatory'),] > 1, 2, any))), 'Mito-hi')
ri = res[c1 %in% cells.use][c2 %in% cells.use][epadj < .05][group1 != group2]
ri = ri[,.SD[which.min(epadj)],lig]
ri = ri[,.SD[which.max(ratio)],.(c1,c2)]
pi = nice_network(edges=ri[,.(source=c1, target=c2, weight=1, label=paste(lig, rec, sep='\n'))], node_colors=f2g, layout='mds')
save_plot(pi, file='Fig_6D.ligrec_network.pdf', nrow=2, ncol=2)


# Table S - Write receptor-ligand interactions
# --------------------------------------------
res_out = res[,.(ident1=c1, ident2=c2, lineage1=group1, lineage2=group2, gene1=gene1, gene2=gene2, mean1=g1, mean2=g2, stat=true, null=shuf, pval, epval, epadj, "lig"=(gene1==lig))]
write.table(res_out, file='../table/Table_S.ligrec_cci.tsv', sep='\t', quote=F, row.names=F)


# Figure 6E - representative images for follicle/inflammatory hits
# ----------------------------------------------------------------

# define cell-cell interactions to plot
ri.use = list(c('Follicular helper T cells', 'B cells', 'CXCL13', 'CXCR5'),
              c('Cycling T cells', 'CD63+CD81+ macrophage', 'IFNG', 'IFNGR2'),
	      c('CD63+CD81+ macrophage', 'Inflammatory fibroblasts', 'IL1B', 'IL1R1'))
ri.use = do.call(rbind, sapply(ri.use, function(a) res[c1 == a[[1]]][c2 == a[[2]]][lig == a[[3]]][rec == a[[4]]], simplify=F))
sa.use = c('V10S15-054_D', 'V10A14-143_C', 'V10A14-143_D')

# plot representative examples
ps = sapply(1:nrow(ri.use), function(i) {
    c1 = ri.use[i]$c1
    c2 = ri.use[i]$c2
    g1 = ri.use[i]$gene1
    g2 = ri.use[i]$gene2
    si = sa.use[[i]]
    scores = data.frame(p[,c1]*p[,c2], vis$data[g1,]*vis$data[g2,])
    colnames(scores) = c(paste(c1, c2, sep='\n'), paste(g1, g2, sep='\n'))
    plot_tsne(vis, scores=scores, ident=F, ymax=.999, pt.size=.75, s.use=si, num_col=1, do.sort=T, cblind=CB)
}, simplify=F)

ps = plot_grid(plotlist=ps, nrow=1)
save_plot(ps, file='Fig_6E.ligrec_vis.pdf', nrow=2, ncol=2)


# ----------------------------------------------------
# FIGURE 6FG - cell-cell interactions with lasso model
# ----------------------------------------------------

# Load lasso model results
# ------------------------
# bi = lasso coefficient
# ri = lasso residuals
# vi = lasso explained variance
bi = readRDS('../3.Lasso/beta.rds')
ri = readRDS('../3.Lasso/residuals.rds')
vi = 1 - readRDS('../3.Lasso/variance.rds')

# Filter hits by explained variance (min 1% varexp)
vmin = 0.01
d = as.data.table(t(apply(which(vi > vmin, arr.ind=T), 1, function(a) c(strsplit(rownames(vi)[[a[[1]]]],'\\.')[[1]], colnames(vi)[[a[[2]]]]))))
colnames(d) = c('lcell', 'rcell', 'gene')
d$lgroup = f2g[d$lcell]
d$varexp = vi[vi > vmin]
d$varmin = apply(vi, 2, min)[d$gene]

# Merge cell-cell interaction data
# --------------------------------

# add cell-cell interaction
D = as.data.table(as.data.frame(bi) %>% rownames_to_column('ct') %>% gather(gene, beta, -ct))[beta != 0][, lcell := gsub('\\..*', '', ct)][, rcell := gsub('.*\\.', '', ct)]
d = merge(d, D, by=c('lcell', 'rcell', 'gene'))
sd = apply(bi, 2, sd)
d[, beta.sd := sd[gene]]

# add expression information
genes.use = sort(unique(d$gene))
D = as.data.table(as.data.frame(as.matrix(gep[, genes.use])) %>% rownames_to_column('lcell') %>% gather(gene, mean_exp, -lcell))
d = merge(d, D, by=c('lcell', 'gene'))
D = apply(gep, 2, max)
d[, max_exp := D[gene]]


# Figure 6F - inflammatory fibroblast cell-cell interactions
# ----------------------------------------------------------
u = d[grepl('Inflammatory', lcell)][grep('Fib|Glia|Peri', rcell, invert=T)][, beta := beta - 1][beta > 0.5][mean_exp > .5]
v = data.frame(u[,.(rcell, gene, beta)] %>% spread(gene, beta), row.names=1)
p1 = ggheatmap(t(v), replace_na=0, vmax=10, legend.title='Beta')
save_plot(p1, file='Fig_6F.iaf_lasso_heatmap.pdf', nrow=1.25, ncol=1)


# Figure 6G - representative examples of IAF interactions
# -------------------------------------------------------

znkg = readRDS('../1.BayesPrism/total.Znkg.rds')

plot_ix = function(di, sa.use=NULL, sa.rmv=NULL, do.log=FALSE, vmax=NA, ymax=.999, ...){

    # get interaction data
    c1 = di$lcell[[1]]
    c2 = di$rcell[[1]]
    gi = di$gene
    
    # average gene expression
    f = function(x){
        x = sign(x)*log2(abs(x)+1)
	2**rowMeans(x)-1
    }
    
    # build dataframe for plot
    qi = data.frame(obs=f(y[,gi,drop=F]), exp=f(znkg[,gi,drop=F]), res=f(y[,gi,drop=F] - znkg[,gi,drop=F]), c1=p[,c1], c2=q[,c2], ix=p[,c1]*q[,c2])
    colnames(qi) = c('Observed', 'Predicted', 'Residual', c1, c2, 'Interaction zone')
    
    # fix column names
    if(do.log){qi[,1:3] = sign(qi[,1:3])*log2(abs(qi[,1:3])+1)}
    
    # plot cell-cell interactions
    plot_tsne(vis, scores=qi, do.sort=T, ident=F, s.use=sa.use, theme_void=T, ymin=1-ymax, ymax=ymax, num_col=3, cblind=CB)
}

# plot representative examples of IAF interactions
u = d[grepl('Inflammatory', lcell)][grep('Fib|Glia|Peri', rcell, invert=T)][, beta := beta - 1][beta > 0.5][mean_exp > .5]
u = u[,.SD[which.max(beta)],gene]
p2 = plot_ix(u[rcell == 'Epithelial'], sa.use='V10A14-143_D', ymax=.999, do.log=F)
p3 = plot_ix(u[rcell == 'T cell'], sa.use='V11Y24-011_C', ymax=.999, do.log=F)
ps = plot_grid(p2,p3,nrow=2)
save_plot(ps, file='Fig_6G.iaf_lasso_residual.pdf', nrow=2, ncol=1.25)


# -----------------------------------------------------------
# Figure 6H - zoom in of inflammatory fibroblast interactions
# -----------------------------------------------------------

# find examples of inflammatory fibroblast IL11 and IL24 interactions
u1 = sort(p[,'Inflammatory fibroblasts']*q[,'Epithelial']*vis$data['IL24',], dec=T)
u2 = sort(p[,'Inflammatory fibroblasts']*q[,'T cell']*vis$data['IL11',], dec=T)

# get spatial coordinates of interacting cell types
c1 = 'V10A14-143_D.TCCCTTGTCTGAAACT'
c2 = 'V11Y24-011_C.TGCGGAGTAAAGGTGC'
s1 = gsub('\\..*', '', c1)
s2 = gsub('\\..*', '', c2)
c1 = as.integer(vis$impos[[s1]][c1,])
c2 = as.integer(vis$impos[[s2]][c2,])

# define zoom coordinates
k = 250
p1 = plot_tsne(vis, 'IL24', scores=data.frame('Inflammatory fibroblast'=p[,'Inflammatory fibroblasts'], 'Epithelial'=q[,'Epithelial']), image=vis$image[[s1]], coords=vis$impos[[s1]], ident=F, xzoom=c(c1[[1]]-k, c1[[1]]+k), yzoom=c(c1[[2]]-k, c1[[2]]+k), pt.size=1.5, num_col=3, cblind=CB)
p2 = plot_tsne(vis, 'IL11', scores=data.frame('Inflammatory fibroblast'=p[,'Inflammatory fibroblasts'], 'Myeloid'=q[,'T cell']), image=vis$image[[s2]], coords=vis$impos[[s2]], ident=F, xzoom=c(c2[[1]]-k, c2[[1]]+k), yzoom=c(c2[[2]]-k, c2[[2]]+k), pt.size=1.5, num_col=3, cblind=CB)
p3 = plot_tsne(vis, 'PDPN', scores=data.frame('Inflammatory fibroblast'=p[,'Inflammatory fibroblasts'], 'Epithelial'=q[,'Epithelial']), image=vis$image[[s1]], coords=vis$impos[[s1]], ident=F, xzoom=c(c1[[1]]-k, c1[[1]]+k), yzoom=c(c1[[2]]-k, c1[[2]]+k), pt.size=1.5, num_col=3, cblind=CB)
ps = plot_grid(p1,p2,p3,nrow=3,align='hv')
save_plot(ps, file='Fig_S5E.iaf_interactions.pdf', nrow=2.25, ncol=2)


# -----------
# Figure 7EF: 
# -----------


# Figure 7E: network plot of interactions between risk genes
# ----------------------------------------------------------
u = d[gene %in% c(gwas.ibd)][,.SD[which.max(varexp)],gene]
v = u[,.(source=lcell, target=rcell, label=gene, weight=1)]
v = v[,list(label = paste(label, collapse='\n'), weight = 1),.(source, target)]
node_colors = sapply(unstack(v[,.(Source=source, Context=target)] %>% gather()), unique)
ps = nice_network(edges=v, node_colors=node_colors)
save_plot(ps, file='Fig_7E.GWAS_lasso_network.pdf', nrow=1.5, ncol=1.5)


# Figure 7F: plot epithelial-tcell interactions
# ---------------------------------------------
w = u[lcell == 'Enterocyte-ANPEP-hi'][rcell == 'T cell']
p1 = plot_ix(w, sa.use='V10S15-055_B', ymax=.999, do.log=F)
save_plot(p1, file='Fig_7F.GWAS_example.pdf', nrow=1.25, ncol=1.5)


# ----------------------------------------------
# FIGURE 7A-D: spatial mapping of IBD risk genes
# ----------------------------------------------


# Figure 7A: GWAS heatmap
# -----------------------

# map ibd risk genes to single-cell and spatial data
# --------------------------------------------------

# map risk genes to highest expressing cell types (single-cell data)
gwas2ident = gwas2group = setNames(rep(NA, length(gwas.ibd)), gwas.ibd)
i = gwas.ibd %in% colnames(gep)
gwas2ident[i] = apply(gep[, gwas.ibd[i]], 2, function(a) names(which.max(a)))
gwas2group[i] = f2g[gwas2ident[i]]

# map risk genes to spatial clusters (spatial data)
scores.cd = score_cells(vis, c('TNF', 'TNFRSF1A', 'TNFRSF1B', gwas.cd), groups=vis$meta.data$cluster)


# plot single-cell and spatial maps
# ---------------------------------

# cluster risk genes by spatial expression
hh = hclust(dist(t(scale(scores.cd))))
Colv = colnames(scores.cd)[hh$order]
col = setNames(ifelse(Colv %in% names(gwas2group), gwas2group[Colv], 'Other'), Colv)

# split risk genes into spatial clusters
gwas.k = cutree(hh, k=6)
gwas.k = tapply(names(gwas.k), gwas.k, c)
names(gwas.k) = c('Follicle', 'Epithelial', 'Inflammatory', 'Muscle', 'Stromal', 'Neuron')
scores.k = score_cells(vis, gwas.k, combine_genes='scale2')

# plot heatmap with dendrogram and cell groups
di = dendro_data_k(hh, 6)
p0 = plot_ggdendro(di, direction='tb', scale.color=set.colors, branch.size=.5, do.label=F) + theme_void() + scale_x_continuous(limits=c(0,length(col)+1), expand=c(0,0), labels=Colv, breaks=1:length(Colv))
p1 = ggheatmap(scale(scores.cd[,Colv]), Colv='none') + scale_x_continuous(limits=c(0,length(col)+1), expand=c(0,0), labels=Colv, breaks=1:length(Colv))
p2 = simple_scatter(1:length(Colv), rep(1, length(Colv)), col=col, pal=set.colors, size=2) + scale_x_continuous(limits=c(0,length(col)+1), expand=c(0,0), labels=Colv, breaks=1:length(Colv)) + theme_void()
l1 = get_legend(p1)
l2 = get_legend(p2)
ps1 = plot_grid(p0, p1+theme(legend.position='none'), p2+theme(legend.position='none'), rel_heights=c(.25, 1, .25), align='v', nrow=3)
ps2 = plot_grid(l1, l2, nrow=2)

ps = plot_grid(ps1, ps2, rel_widths=c(1,.25), nrow=1)
save_plot(ps, file='Fig_7A.GWAS_heatmap.pdf', nrow=1.25, ncol=2)


# Figure 7B: risk gene enrichment in visium clusters
# --------------------------------------------------

# A) risk gene enrichment in spot clusters

# scale the log-proportions of cell groups for regression
qq = scale(log(q))

if(FALSE){

    # generate "background" gene sets for null distribution
    gwas.bg = sapply(gwas.k, function(a) select_bg_genes(vis, genes.use=a, nbins=25, nperm=1000), simplify=F)
    
    # score expression of risk genes and background gene sets
    scores.k = score_cells(vis, gwas.k, groups=vis$meta.data$cluster, combine_genes='scale2')
    scores.bg = sapply(1:1000, function(i) score_cells(vis, sapply(gwas.bg, function(a) a[,i]), groups=vis$meta.data$cluster, combine_genes='scale2'), simplify=F)
    
    # estimate empirical p-values and save results
    pvals = Reduce('+', sapply(scores.bg, function(a) a >= scores.k, simplify=F))

    # adjust p-values
    pvals = (pvals + 1)/1001.
    pvals[] = p.adjust(pvals[], 'fdr')
    
    saveRDS(pvals, file='gwas_enrichment.pvals.rds')
} else {

    # alternatively, load saved results
    pvals = readRDS('gwas_enrichment.pvals.rds')
}

# plot significant enrichments
i = apply(pvals, 1, min) <= .01
p1 = ggheatmap(-log10(pvals[i,]), pal='Blues', vmin=-log10(.05), xlab='CD risk gene clusters', ylab='Visium spot clusters', legend.title='-log10(P)')
save_plot(p1, file='Fig_7B.GWAS_enrichment_heatmap.pdf', nrow=1, ncol=1)


if(FALSE){

    # score expression of risk genes in each spot cluster
    scores.k = score_cells(vis, gwas.k, combine_genes='scale2')
    # estimate residual of regression against cell type proportions
    scores.k = sapply(scores.k, function(a) lm(a ~ qq)$resid)
    # aggregate residuals by taking the mean per spot cluster
    scores.k = nice_agg(scores.k, vis$meta.data$cluster, 'mean')

    # generate "background" gene sets for null distribution
    gwas.bg = sapply(gwas.k, function(a) select_bg_genes(vis, genes.use=a, nbins=25, nperm=1000), simplify=F)
    scores.bg = sapply(1:1000, function(i){
        # score expression of "background" genes in each spot cluster
        scores = score_cells(vis, sapply(gwas.bg, function(a) a[,i]), combine_genes='scale2')
	# estimate residual of regression against cell type proportions
	scores = sapply(scores, function(a) lm(a ~ qq)$resid)
	# aggregate residuals by taking the mean per spot cluster
	scores = nice_agg(scores, vis$meta.data$cluster, 'mean')
	scores
    }, simplify=F)

    # estimate empirical p-values and save results
    pvals = Reduce('+', sapply(scores.bg, function(a) a >= scores.k, simplify=F))
    saveRDS(pvals, file='gwas_enrichment.cell_proportions.pvals.rds')
    
} else {

    # alternatively, load saved results
    pvals = readRDS('gwas_enrichment.cell_proportions.pvals.rds')
}

# adjust p-values and select significant tests
pvals = (pvals+1)/1001
pvals[] = p.adjust(pvals[], 'fdr')
i = apply(pvals, 1, min) <= .05

# plot significant enrichments (controlling for cell type composition)
p1 = ggheatmap(-log10(pvals[i,]), pal='Blues', vmin=-log10(.05), xlab='CD risk gene clusters', ylab='Visium spot clusters', legend.title='-log10(P)')
save_plot(p1, file='Fig_7C.GWAS_enrichment_heatmap2.pdf', nrow=1, ncol=1)


# Figure 7C: risk gene enrichment in spot clusters, controlling for cell type composition
# ---------------------------------------------------------------------------------------

# scale the log-proportions of cell groups for regression
qq = scale(log(q))

# get spot cluster assignments
names.use = c('Epithelial', 'Inflammatory', 'Follicle')
v = vis$meta.data$cluster
levels(v)[levels(v) %in% c('Mature', 'TA', 'Stem')] = 'Epithelial'

# calculate scores and residuals
scores.k = score_cells(vis, gwas.k, combine_genes='scale2')
resids.k = sapply(scores.k, function(a) lm(a ~ qq)$resid)

# boxplot showing risk gene score for each cluster vs. "other"
d = sapply(names.use, function(a) data.frame(type=ifelse(is.na(v) | (v != a), 'Other', a), score=qtrim(scores.k[,a], qmax=.995)), simplify=F)
p1 = sapply(names(d), function(aa){a = d[[aa]]; print(aa)
    ggplot(a, aes(x=type, y=score)) + geom_boxplot(aes(fill=type), outlier.shape=NA) + geom_signif(comparisons=list(c(aa, 'Other'))) + theme_cowplot() + xlab('') + ylab('Mean expression (log2(TP10K+1))') +
    theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1), legend.position='none') + scale_fill_manual('', values=c(set.colors[1:3], '#cccccc'), breaks=c(names.use, 'Other'))
    }, simplify=F)
p1 = plot_grid(plotlist=p1,nrow=1,align='h')

# boxplot showing risk gene residual for each cluster vs. "other"
d = sapply(names.use, function(a) data.frame(type=ifelse(is.na(v) | (v != a), 'Other', a), score=qtrim(resids.k[,a], qmax=.995)), simplify=F)
p2 = sapply(names(d), function(aa){a = d[[aa]]; print(aa)
    ggplot(a, aes(x=type, y=score)) + geom_boxplot(aes(fill=type), outlier.shape=NA) + geom_signif(comparisons=list(c(aa, 'Other'))) + theme_cowplot() + xlab('') + ylab('Residuals (log2(TP10K+1))') +
    theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1), legend.position='none') + scale_fill_manual('', values=c(set.colors[1:3], '#cccccc'), breaks=c(names.use, 'Other'))
    }, simplify=F)
p2 = plot_grid(plotlist=p2,nrow=1,align='h')

# merge plots
ps = plot_grid(p1,p2,nrow=2)
save_plot(ps, file='Fig_7C.GWAS_enrichment_boxplot.pdf', nrow=2, ncol=1.5)


# Figure 7D: plot of spot clusters + risk gene expression
# -------------------------------------------------------

# initialize variables
se = 'V10A14-143_A'
si = 'V10A14-143_D'
sf = 'V11Y24-011_D'

# calculate for each visium spot the mean expression of the risk gene clusters (gwas.k)
scores.k = score_cells(vis, gwas.k, combine_genes='scale2')

# calculate for each visium spot the annotated cell group (collapsing epithelial cells)
v = vis$meta.data$cluster
levels(v)[levels(v) %in% c('Mature', 'TA', 'Stem')] = 'Epithelial'
v = sapply(levels(v), function(a) as.numeric(v == a))
rownames(v) = rownames(vis$meta.data)
v[is.na(v)] = 0

# convert the mean gene expression into rgb (each channel = cell type proportions)
names.use = c('Epithelial', 'Inflammatory', 'Follicle')
w = scores.k[,names.use]
f = apply(w[grep(paste(se, si, sf, sep='|'), rownames(w)),], 2, quantile, ymax)
w = scale(w, center=F, scale=f)

# plot the cell type assignments for each tissue section
p1 = plot_tsne(vis, do.sort=F, pt.size=1, theme_void=T, multichannel=v[,names.use], ident=F, coords=vis$impos[[se]], do.title=F, cblind=CB)
p2 = plot_tsne(vis, do.sort=F, pt.size=1, theme_void=T, multichannel=v[,names.use], ident=F, coords=vis$impos[[si]], do.title=F, cblind=CB)
p3 = plot_tsne(vis, do.sort=F, pt.size=1, theme_void=T, multichannel=v[,names.use], ident=F, coords=vis$impos[[sf]], do.title=F, cblind=CB)

# plot the gene expression scores for each tissue section
p4 = plot_tsne(vis, do.sort=F, pt.size=1, theme_void=T, multichannel=w[,names.use], ident=F, coords=vis$impos[[se]], ymax=.99, do.title=F, multichannel_rescale=F, cblind=CB)
p5 = plot_tsne(vis, do.sort=F, pt.size=1, theme_void=T, multichannel=w[,names.use], ident=F, coords=vis$impos[[si]], ymax=.99, do.title=F, multichannel_rescale=F, cblind=CB)
p6 = plot_tsne(vis, do.sort=F, pt.size=1, theme_void=T, multichannel=w[,names.use], ident=F, coords=vis$impos[[sf]], ymax=.99, do.title=F, multichannel_rescale=F, cblind=CB)

# merge plots (top = visium clusters, bottom = risk genes)
nl = theme(legend.position='none')
lg = get_legend(p1)
ps = lapply(list(p1,p2,p3,p4,p5,p6), function(a) a+nl)
ps = plot_grid(plotlist=ps, nrow=2)
ps = plot_grid(ps,lg,nrow=1,rel_widths=c(1,.15))
save_plot(ps, file='Fig_7D.GWAS_expression_examples.pdf', nrow=1.5, ncol=1.5)
