set.seed(10)
session_info <- sessionInfo(); date_of_run <- Sys.time()

load("../../raw_data/Zeisel_expr.RData") # loads file "dat"
labels <- read.csv("../../raw_data/Zeisel_cell_info.txt", sep = "\t")
dat <- as.matrix(dat)
label_vec <- labels$level1class

##########

dat_count <- dat
label_vec <- as.factor(as.character(label_vec))

# next remove genes
zz <- apply(dat, 2, function(x){length(which(x!=0))})
idx <- which(zz > nrow(dat)/100)
dat <- dat[,idx]
dat_count <- dat_count[,idx]

# normalize
dat <- t(sapply(1:nrow(dat), function(i){10^4 * dat[i,]/sum(dat[i,])}))
rownames(dat) <- 1:nrow(dat)

# try a series of SPCAs
k <- 5
lvls <- 20
v_seq <- exp(seq(log(1), log(sqrt(ncol(dat))/2), length.out = lvls))
res_list <- vector("list", lvls)

spca_func <- function(i){
  res <- PMA::SPC(dat, sumabsv = v_seq[i], K = k, trace = F)
  print(paste0("Finished SPC for level ", i))
  res
}

res_list <- foreach::"%dopar%"(foreach::foreach(i = 1:lvls), spca_func(i))

# run sPCA
spca_mat <- cbind(v_seq, t(sapply(res_list, function(x){c(length(unique(sort(unlist(apply(x$v, 2, function(y){which(y != 0)}))))), x$prop.var.explained[5])})))
idx <- min(which(spca_mat[,2] == ncol(dat)))
target_var <- spca_mat[idx,3]
idx <- min(intersect(which(spca_mat[,2] >= 500), which(spca_mat[,3] >= 0.9*target_var)))
spca_idx <- sort(unlist(apply(res_list[[idx]]$v, 2, function(x){which(x != 0)})))
spca_hvg <- colnames(dat)[spca_idx]

# run DESCEND
res_descend <- descend::runDescend(t(dat_count), n.cores = ncores)
descend_hvg <- descend::findHVG(res_descend, threshold = 50)$HVG.genes

idx <- which(colnames(dat) %in% c(descend_hvg, spca_hvg))
dat <- dat[,idx]

save.image(paste0("../results/step0_zeisel_preprocessing", suffix, ".RData"))

print(paste0(Sys.time(), ": Finished preprocessing"))
source_code_info <- c(source_code_info, readLines("../main_supplement/step0_zeisel_preprocessing.R"))
save.image(paste0("../results/step0_zeisel_preprocessing", suffix, ".RData"))
print(warnings())
