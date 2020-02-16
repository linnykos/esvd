set.seed(10)
load(paste0("../results/step2_naive_svd", suffix, ".RData"))

max_val <- 5000
scalar_vec <- 1:5
res_list <- vector("list", length(scalar_vec))

for(i in 1:length(scalar_vec)){
  print(paste0("Trying scalar value = ", scalar_vec[i]))
  init <- eSVD::initialization(dat_impute_NA, family = family,
                                       k = k, max_val = max_val)
  res_list[[i]] <- eSVD::fit_factorization(dat_impute_NA, u_mat = init$u_mat, v_mat = init$v_mat,
                                                   family = family, reparameterize = T,
                                                   max_iter = 25, max_val = max_val,
                                                   scalar = scalar_vec[i],
                                                   return_path = F, cores = ncores,
                                                   verbose = T)
  save.image(paste0("../results/step3_scalar_heuristic", suffix, "_tmp.RData"))
}

.l2norm <- function(x){sqrt(sum(x^2))}
quality_vec <- sapply(res_list, function(x){
  pred_mat <- 1/(x$u_mat %*% t(x$v_mat))

  mat <- cbind(dat_impute[missing_idx], pred_mat[missing_idx])

  pca_res <- stats::princomp(mat)
  diag_vec <- c(1,1)

  # plot(mat[,1], mat[,2], asp = T, pch = 16, col = rgb(0,0,0,0.2)); lines(c(0,1e6), c(0,1e6), col = "red", lwd = 2); lines(c(0, 1e6), c(0, 1e6*pca_res$loadings[2,1]/pca_res$loadings[1,1]), col = "blue", lwd = 2, lty = 2)

  acos(diag_vec %*% pca_res$loadings[,1] / (.l2norm(diag_vec) * .l2norm(pca_res$loadings[,1])))
})

scalar_val <- scalar_vec[which.min(quality_vec)]

rm(list = c("i", "init"))
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




