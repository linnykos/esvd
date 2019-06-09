library(descend)
load("../../raw_data/marques.RData")

dat <- marques$counts

cell_types <- unique(marques$cell.info$cell.type)

set.seed(10)
cell_idx <- unlist(lapply(cell_types, function(x){
  tmp <- which(marques$cell.info$cell.type == x)
  sample(tmp, round(length(tmp)/2.5))
}))
dat <- dat[cell_idx,]
dim(dat)

# # remove genes with too many 0's
zz <- apply(dat, 2, function(x){length(which(x!=0))})
dat <- dat[,which(zz > 30)]

# try a series of SPCAs
k <- 5
lvls <- 20
v_seq <- exp(seq(log(1), log(log(ncol(dat))), length.out = lvls))
res_list <- vector("list", lvls)

spca_func <- function(i){
  res <- PMA::SPC(dat, sumabsv = v_seq[i], K = k, trace = F)
  print(paste0("Finished SPC for level ", i))
  res
}

doMC::registerDoMC(cores = ncores)
res_list <- foreach::"%dopar%"(foreach::foreach(i = 1:lvls), spca_func(i))

print(paste0(Sys.time(), ": Finished SPC"))

# run DESCEND
res_descend <- descend::runDescend(t(dat), n.cores = ncores)

rm(list = c("idx", "zz", "k", "lvls", "reorder_idx", "column_vec"))
print(paste0(Sys.time(), ": Finished screening"))
save.image(paste0("../results/step0_screening", suffix, ".RData"))
