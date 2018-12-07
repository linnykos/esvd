set.seed(10)
load("../results/step2_naive_svd.RData")

max_val <- 5000
scalar_vec <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 100)
res_list <- vector("list", length(scalar_vec))

for(i in 1:length(scalar_vec)){
  init <- singlecell::initialization(dat_impute_NA, family = family, extra_weight = extra_weight,
                                       k = k, max_val = max_val)
  res_list[[i]] <- singlecell::fit_factorization(dat_impute_NA, u_mat = init$u_mat, v_mat = init$v_mat,
                                                   family = family, reparameterize = T,
                                                   max_iter = 25, max_val = max_val,
                                                   scalar = scalar_vec[i], extra_weight = extra_weight,
                                                   return_path = F, cores = 15,
                                                   verbose = T)
  save.image(paste0("../results/step3_scalar_heuristic", suffix, "_tmp.RData"))
}

.l2norm <- function(x){sqrt(sum(x^2))}
quality_vec <- sapply(res_list, function(x){
  pred_mat <- 1/(x$u_mat %*% t(x$v_mat))
  pred_mat <- t(sapply(1:nrow(pred_mat), function(x){
    pred_mat[x,] * extra_weight[x]
  }))

  mat <- cbind(dat_impute[missing_idx], pred_mat[missing_idx])
  mat <- mat[which(mat[,1] <= 1800),]

  pca_res <- stats::princomp(mat)
  diag_vec <- c(1,1); diag_vec <- diag_vec/.l2norm(diag_vec)

  # plot(mat[,1], mat[,2], asp = T, pch = 16, col = rgb(0,0,0,0.2)); lines(c(0,1e6), c(0,1e6), col = "red", lwd = 2); lines(c(0, 1e6), c(0, 1e6*pca_res$loadings[2,1]/pca_res$loadings[1,1]), col = "blue", lwd = 2, lty = 2)

  acos(diag_vec %*% pca_res$loadings[,1])
})

scalar_val <- scalar_vec[which.min(quality_vec)]

rm(list = c("i", "init"))
print(paste0(Sys.time(), ": Finished scalar heuristic"))
save.image(paste0("../results/step3_scalar_heuristic", suffix, ".RData"))


# for(i in 1:length(res_list)){
#   plot(4/pred_mat[idx], dat_impute[idx], asp = T, main = i, pch = 16,
#        col = rgb(0,0,0,0.2))
#   lines(c(-1e10,1e10), c(-1e10,1e10), col = "red", lwd = 2)
# }
# ############
#
# plot(pred_naive[idx], dat_impute[idx], asp = T, main = i, pch = 16,
#      col = rgb(0,0,0,0.2))
# lines(c(-1e10,1e10), c(-1e10,1e10), col = "red", lwd = 2)


