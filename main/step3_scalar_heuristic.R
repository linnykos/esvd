set.seed(10)
load("../results/step2_naive_svd.RData")

# max_val <- 2000
max_val <- 8
scalar_vec <- c(0.5, 1, 1.5, 1.75, 2, 2.25, 2.5, 3, 5)
res_list <- vector("list", length(scalar_vec))
extra_weight <- apply(dat_impute, 1, mean)

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

quality_vec <- sapply(res_list, function(x){
  pred_mat <- x$u_mat %*% t(x$v_mat)
  pred_mat <- t(sapply(1:nrow(pred_mat), function(x){
    pred_mat[x,] * extra_weight[x]
  }))
  pred_mat <- exp(pred_mat)

  mat <- cbind(dat_impute[missing_idx], pred_mat[missing_idx])

  # plot(mat[,1], mat[,2], asp = T, pch = 16, col = rgb(0,0,0,0.2)); lines(c(0,1e6), c(0,1e6), col = "red", lwd = 2)

  pca_res <- stats::princomp(mat)
  diag_vec <- c(1,1); diag_vec <- diag_vec/.l2norm(diag_vec)

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


