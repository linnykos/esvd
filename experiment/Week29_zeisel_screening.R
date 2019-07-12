rm(list=ls())
load("../../SOUPR/data/zeisel.rda")

dat <- zeisel$counts

idx <- grep("oligodendrocytes", zeisel$cell.info$cell.type)
dat <- dat[idx,]
dim(dat)

# first select some genes based on SPCA
cov_dat <- stats::cov(dat)
