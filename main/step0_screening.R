rm(list=ls())
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
lvls <- 10
v_seq <- exp(seq(log(1), log(sqrt(ncol(dat))), length.out = lvls))
res_list <- vector("list", lvls)

for(i in 1:lvls){
  print(i)
  save.image(file = "../results/step0_screening_tmp.RData")
  res_list[[i]] <- PMA::SPC(dat, sumabsv = v_seq[i], K = k, trace = F)
}

save.image("../results/step0_screening.RData")
