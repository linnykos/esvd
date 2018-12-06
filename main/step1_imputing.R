set.seed(10)
load("../results/step0_screening.RData")

extra_weight <- apply(dat, 1, mean)

res_hvg <- descend::findHVG(res_descend, threshold = 12)
length(res_hvg$HVG.genes)

idx1 <- sort(unlist(apply(res_list[[5]]$v, 2, function(x){which(x != 0)})))
idx2 <- which(colnames(dat) %in% res_hvg$HVG.genes)
idx <- sort(unique(c(idx1, idx2)))

dat <- dat[,idx]
dim(dat)

dropout_mat <- singlecell::dropout(dat)
zero_mat <- singlecell::find_true_zeros(dropout_mat, num_neighbors = 50)
idx <- which(is.na(zero_mat))

dat_impute <- singlecell::scImpute(dat, drop_idx = idx, Kcluster = 5,
                                     verbose = F, weight = 1)

rm(list = c("idx1", "idx2", "idx"))
print(paste0(Sys.time(), ": Finished imputing"))
save.image("../results/step1_imputing.RData")
