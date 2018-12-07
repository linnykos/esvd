set.seed(10)
load("../results/step0_screening.RData")

reweight_factor <- rowSums(dat)
extra_weight <- rep(1, nrow(dat))

res_hvg <- descend::findHVG(res_descend, threshold = 12)
length(res_hvg$HVG.genes)

idx1 <- sort(unlist(apply(res_list[[5]]$v, 2, function(x){which(x != 0)})))
idx2 <- which(colnames(dat) %in% res_hvg$HVG.genes)
idx <- sort(unique(c(idx1, idx2)))

dat <- dat[,idx]
dat <- t(sapply(1:nrow(dat), function(i){dat[i,]/reweight_factor[i]}))
dat <- log(dat+1)
dat <- dat * 10/mean(dat)
dim(dat)

dropout_mat <- singlecell::dropout(dat)
zero_mat <- singlecell::find_true_zeros(dropout_mat, num_neighbors = 200)
idx <- which(is.na(zero_mat))

dat_impute <- singlecell::scImpute(dat, drop_idx = idx, Kcluster = 5,
                                     verbose = T, weight = 1)

# png("../figure/main/data.png", height = 2400, width = 1000, res = 300, units = "px")
# par(mar = rep(0.5, 4))
# singlecell:::.plot_singlecell(dat)
# graphics.off()
#
# png("../figure/main/data_impute.png", height = 2400, width = 1000, res = 300, units = "px")
# par(mar = rep(0.5, 4))
# singlecell:::.plot_singlecell(dat_impute)
# graphics.off()

rm(list = c("idx1", "idx2", "idx", "res_descend", "res_list", "v_seq", "k"))
print(paste0(Sys.time(), ": Finished imputing"))
save.image("../results/step1_imputing.RData")
