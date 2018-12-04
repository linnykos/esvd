set.seed(10)
load("../results/step2_naive_svd.RData")

k <- 5
max_val <- 2000
scalar_vec <- c(1, 1.5, 1.75, 2, 2.25, 2.5, 3)
res_list <- vector("list", length(scalar_vec))

for(i in 1:length(scalar_vec)){
  init <- singlecell::initialization(dat_impute_NA, family = "gaussian", scalar = scalar_vec[i],
                                       k = k, max_val = max_val)
  res_list[[i]] <- singlecell::fit_factorization(dat_impute_NA, u_mat = init$u_mat, v_mat = init$v_mat,
                                                   family = "gaussian",  reparameterize = T,
                                                   max_iter = 25, max_val = max_val,
                                                   scalar = scalar_vec[i],
                                                   return_path = F, cores = 15,
                                                   verbose = T)
  save.image("../results/step2b_scalar_heuristic_tmp.RData")
}

quality_vec <- sapply(res_list, function(x){
  pred_mat <- 4/(x$u_mat %*% t(x$v_mat))
  mat <- cbind(dat_impute[idx], pred_mat[idx])

  pca_res <- stats::princomp(mat)
  diag_vec <- c(1,1); diag_vec <- diag_vec/.l2norm(diag_vec)

  acos(diag_vec %*% pca_res$loadings[,1])
})

scalar_val <- scalar_vec[which.min(quality_vec)]

rm(list = c("i", "init"))
print(paste0(Sys.time(), ": Finished scalar heuristic"))
save.image("../results/step3_scalar_heuristic.RData")


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


