rm(list=ls())
library(descend)
load("../../SOUPR/data/zeisel.rda")

dat <- zeisel$counts

idx <- grep("oligodendrocytes", zeisel$cell.info$cell.type)
dat <- dat[idx,]
dim(dat)

# remove genes with too many 0's
zz <- apply(dat, 2, function(x){length(which(x == 0))/length(x)})
idx <- which(zz <= .95)
dat <- dat[,idx]
dim(dat)

# try a series of SPCAs
k <- 5
lvls <- 20
v_seq <- exp(seq(log(1), log(log(ncol(dat))), length.out = lvls))
res_list <- vector("list", lvls)

for(i in 1:lvls){
  print(paste0(Sys.time(), ": Level ", i))
  save.image(file = "../results/step0_screening_tmp.RData")
  res_list[[i]] <- PMA::SPC(dat, sumabsv = v_seq[i], K = k, trace = F)
}

# run DESCEND
res_descend <- descend::runDescend(dat, n.cores = 10)
res_hvg <- descend::findHVG(res_descend)

save.image("../results/step0_screening.RData")
