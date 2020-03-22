set.seed(10)
load(paste0("../results/step0_screening", suffix, ".RData"))

res_hvg <- descend::findHVG(res_descend, threshold = 50)
length(res_hvg$HVG.genes)

spca_summary <- cbind(v_seq, t(sapply(res_list, function(x){c(length(unique(sort(unlist(apply(x$v, 2, function(y){which(y != 0)}))))), x$prop.var.explained[5])})))
idx <- min(which(spca_mat[,2] == max(spca_mat[,2])))
target_var <- spca_mat[idx,3]
idx <- min(intersect(which(spca_mat[,2] >= 500), which(spca_mat[,3] >= 0.9*target_var)))
spca_idx <- sort(unlist(apply(res_list[[idx]]$v, 2, function(x){which(x != 0)})))

descend_idx <- which(colnames(dat) %in% res_hvg$HVG.genes)
gene_idx <- sort(unique(c(spca_idx, descend_idx)))
dat <- dat[,gene_idx]

dat_impute <- dat

reweight_factor <- rowSums(dat_impute)
dat_impute <- t(sapply(1:nrow(dat_impute), function(i){dat_impute[i,]/reweight_factor[i]}))
dat_impute <- dat_impute * 1000/max(dat_impute)
dim(dat_impute)

###########

# png("../figure/main/data_impute.png", height = 2400, width = 1000, res = 300, units = "px")
# par(mar = rep(0.5, 4))
# eSVD:::.plot_singlecell(dat_impute)
# graphics.off()

rm(list = c("spca_idx", "descend_idx", "gene_idx", "res_descend",
            "res_hvg", "res_list", "k", "tmp", "dat", "reweight_factor", "idx"))
print(paste0(Sys.time(), ": Finished imputing"))
save.image(paste0("../results/step1_imputing", suffix, ".RData"))
