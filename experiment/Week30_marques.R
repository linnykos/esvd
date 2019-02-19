## svd of imputed marques
rm(list=ls())
load("../results/tmp.RData")
svd_dat <- zz

yy <- svd_dat$u[,1:2] %*% diag(sqrt(svd_dat$d[1:2]))

load("../../raw_data/marques.RData")

idx <- grep("MO", marques$cell.info$cell.type)
tab <-  marques$cell.info[idx,]
col_idx <- as.numeric(as.factor(as.character(tab[,2])))

png("../figure/experiment/Week30_marques_MO.png", height = 1200, width = 1200,
    res = 300, units = "px")
plot(yy[,1], yy[,2], col = col_idx, pch = 16)
graphics.off()
