rm(list=ls())
load("../../SOUP/data/zeisel.rda")

dat <- zeisel$counts

dat <- t(apply(dat, 1, function(x){x/sum(x)}))
dat <- log10(dat + 1)

idx <- which(colnames(dat) %in% zeisel$select.genes)
dat <- dat[,idx]
dim(dat)

table(zeisel$cell.info$cell.type)
idx <- grep("oligodendrocytes", zeisel$cell.info$cell.type)
dat_subset <- dat[idx,]
dim(dat_subset)

x <- dat_subset[,17]
zz <- .em_mixture(x)
.hist_augment(x, param_list = list(zz), multiplier = 0.02, lwd = 3)
zz <- .em_mixture(x, mixture = "exponential.tgaussian")
.hist_augment(x, param_list = list(zz), multiplier = 0.02, lwd = 3)

x <- dat_subset[,340]
zz <- .em_mixture(x)
.hist_augment(x, param_list = list(zz), multiplier = 1, lwd = 3)
zz <- .em_mixture(x, mixture = "exponential.tgaussian")
.hist_augment(x, param_list = list(zz), multiplier = 1, lwd = 3)

x <- dat_subset[,1900]
zz <- .em_mixture(x)
.hist_augment(x, param_list = list(zz), multiplier = 0.1, lwd = 3)
zz <- .em_mixture(x, mixture = "exponential.tgaussian")
.hist_augment(x, param_list = list(zz), multiplier = 0.1, lwd = 3)

