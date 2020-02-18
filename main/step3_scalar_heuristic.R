set.seed(10)
load(paste0("../results/step2_naive_svd", suffix, ".RData"))

max_val <- 5000
scalar_vec <- seq(1, 5, by = 0.5)
res_list <- vector("list", length(scalar_vec))
cv_trials <- 10
n <- nrow(dat_impute); d <- ncol(dat_impute)

for(i in 1:length(scalar_vec)){
  print(paste0("Trying scalar value = ", scalar_vec[i]))
  res_list[[i]] <- vector("list", cv_trials)

  for(j in 1:cv_trials){
    print(paste0("On trial ", j))
    # set missing values
    set.seed(10*j)
    missing_idx <- rbind(do.call(rbind, (lapply(1:n, function(x){
      cbind(x, sample(1:d, 4))
    }))), do.call(rbind, (lapply(1:d, function(x){
      cbind(sample(1:n, 4), d)
    }))))

    dat_impute_NA <- dat_impute
    for(i in 1:nrow(missing_idx)){
      dat_impute_NA[missing_idx[i,1], missing_idx[i,2]] <- NA
    }
    missing_idx <- which(is.na(dat_impute_NA))

    # fit
    init <- eSVD::initialization(dat_impute_NA, family = family,
                                 k = k, max_val = max_val)
    res_list[[i]][[j]] <- eSVD::fit_factorization(dat_impute_NA, u_mat = init$u_mat, v_mat = init$v_mat,
                                             family = family,
                                             max_iter = 100, max_val = max_val,
                                             scalar = scalar_vec[i],
                                             return_path = F, cores = ncores,
                                             verbose = T)
    save.image(paste0("../results/step3_scalar_heuristic", suffix, "_tmp.RData"))
  }
}

.l2norm <- function(x){sqrt(sum(x^2))}
quality_vec <- sapply(res_list, function(x){

  mat_list <- lapply(x, function(y){
    pred_mat <- 1/(y$u_mat %*% t(y$v_mat))

    cbind(dat_impute[missing_idx], pred_mat[missing_idx])
  })

  mat <- do.call(rbind, mat_list)

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




