rm(list=ls())
load("../../SOUP/data/zeisel.rda")

dat <- zeisel$counts

dat <- t(apply(dat, 1, function(x){x/sum(x)}))
dat <- log10(dat + 1.01)

idx <- which(colnames(dat) %in% zeisel$select.genes)
dat <- dat[,idx]
dim(dat)

table(zeisel$cell.info$cell.type)
idx <- grep("oligodendrocytes", zeisel$cell.info$cell.type)
dat_subset <- dat[idx,]
dim(dat_subset)

source("../experiment/em_gamma_normal.R")
# diff_vec <- sapply(1:ncol(dat_subset), function(i){
#   if(i %% floor(ncol(dat_subset)/10) == 0) cat('*')
#   x <- dat_subset[,i]
#   tryCatch({
#     res <- get_mix(x, prop_init = 0.01)
#     res[4] - res[2]/res[3]
#   }, error = function(e) {
#     NA
#   }
#   )
# })

diff_vec <- sapply(1:ncol(dat_subset), function(i){
  if(i %% floor(ncol(dat_subset)/10) == 0) cat('*')
  x <- dat_subset[,i]
  zz <- hist(x, breaks = 100, plot = F)
  sum(zz$counts[-1])/sum(zz$counts)
})

quantile(diff_vec, na.rm = T)
range(diff_vec, na.rm = T)
