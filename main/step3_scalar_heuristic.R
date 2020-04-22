set.seed(10)
load(paste0("../results/step2_naive_svd", suffix, ".RData"))

paramMat_esvd <- as.matrix(expand.grid(c(0.5, 1, 2, 4), c(3,5,10)))
# paramMat_esvd <- as.matrix(expand.grid(c(500, 1000, 5000, 10000, 50000), c(5,10,20,30)))
colnames(paramMat_esvd) <- c("scalar", "k")
esvd_missing_list <- vector("list", nrow(paramMat_esvd))
fitting_distr <- "curved_gaussian"

for(i in 1:nrow(paramMat_esvd)){
  print(paste0("On parameter setting row ", i))
  tmp_list <- vector("list", length(cv_trials))

  for(j in 1:cv_trials){
    print(paste0("On trial ", j))

    # set missing values
    dat_impute_NA <- dat_impute
    dat_impute_NA[missing_idx_list[[j]]] <- NA

    # fit
    set.seed(10)
    init <- eSVD::initialization(dat_impute_NA, family = fitting_distr,
                                 k = paramMat_esvd[i,"k"], max_val = max_val, scalar = paramMat_esvd[i,"scalar"])
    tmp_list[[j]] <- eSVD::fit_factorization(dat_impute_NA, u_mat = init$u_mat, v_mat = init$v_mat,
                                   family = fitting_distr,
                                   max_iter = 50, max_val = max_val,
                                   scalar = paramMat_esvd[i,"scalar"],
                                   return_path = F, cores = ncores,
                                   verbose = T)
    save.image(paste0("../results/step3_scalar_heuristic", suffix, "_tmp.RData"))
  }

  esvd_missing_list[[i]] <- tmp_list
}


rm(list = c("j", "i", "init", "tmp_list", "dat_impute_NA"))
print(paste0(Sys.time(), ": Finished scalar heuristic"))
save.image(paste0("../results/step3_scalar_heuristic", suffix, ".RData"))


# idx <- which.min(quality_vec)
# pred_mat <- 1/(res_list[[idx]]$u_mat %*% t(res_list[[idx]]$v_mat))
# tmp <- cbind(pred_mat[missing_idx], dat_impute[missing_idx])
# pca_res <- princomp(tmp)
# plot(pred_mat[missing_idx], dat_impute[missing_idx], asp = T, pch = 16,
#      col = rgb(0,0,0,0.2), main = "Our embedding",
#      xlab = "Predicted value", ylab = "Observed and masked value")
# lines(c(-1e10,1e10), c(-1e10,1e10), col = "red", lwd = 2)
# lines(c(-1e10,1e10)*pca_res$loadings[1,1], c(-1e10,1e10)*pca_res$loadings[1,2],
#       col = "blue", lwd = 2, lty = 2)
#
#
# ##############
#
# tmp <- cbind(pred_naive[missing_idx], dat_impute[missing_idx])
# pca_res <- princomp(tmp)
# plot(pred_naive[missing_idx], dat_impute[missing_idx], asp = T, pch = 16,
#      col = rgb(0,0,0,0.2), main = "Naive embedding",
#      xlab = "Predicted value", ylab = "Observed and masked value")
# lines(c(-1e10,1e10), c(-1e10,1e10), col = "red", lwd = 2)
# lines(c(-1e10,1e10)*pca_res$loadings[1,1], c(-1e10,1e10)*pca_res$loadings[1,2],
#       col = "blue", lwd = 2, lty = 2)




