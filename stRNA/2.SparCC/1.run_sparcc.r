library(SpiecEasi)
source('../load_data.r')
res = sparcc(p*colSums(vis$counts))
saveRDS(res, 'sparcc.res.rds')
