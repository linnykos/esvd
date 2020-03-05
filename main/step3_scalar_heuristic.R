set.seed(10)
load(paste0("../results/step2_naive_svd", suffix, ".RData"))

# first, determine what the scalar should be
## warning: I'll do a janky tuning procedure for now, we'll fix it to be part of the code later on
set.seed(10)
init <- eSVD::initialization(dat_impute, family = "exponential",
                             k = k, max_val = max_val)
exp_fit <- eSVD::fit_factorization(dat_impute, u_mat = init$u_mat, v_mat = init$v_mat,
                               family = "exponential",
                               max_iter = 50, max_val = max_val,
                               return_path = F, cores = ncores,
                               verbose = T)
scalar_val_init <- eSVD::tuning(dat_impute, exp_fit$u_mat, exp_fit$v_mat, family = "curved_gaussian")

## fit 5 times, each time alternating between scalar_val and fitting
fitting_iter <- 5
scalar_val_vec <- rep(NA, fitting_iter)
scalar_val_vec[1] <- scalar_val_init
for(i in 2:fitting_iter){
  init <- eSVD::initialization(dat_impute, family = "curved_gaussian", k = k, max_val = max_val)
  fit <- eSVD::fit_factorization(dat_impute, u_mat = init$u_mat, v_mat = init$v_mat,
                                     family = "curved_gaussian",  reparameterize = T,
                                     max_iter = 50, max_val = max_val,
                                     scalar = scalar_val_vec[i-1],
                                     return_path = F, cores = ncores,
                                     verbose = T)

  scalar_val_vec[i] <- eSVD::tuning(dat_impute, fit$u_mat, fit$v_mat, family = "curved_gaussian")
}

# next, apply the missing value diagnostic
res_list <- vector("list", length(cv_trials))

for(j in 1:cv_trials){
  print(paste0("On trial ", j))

  # set missing values
  dat_impute_NA <- dat_impute
  dat_impute_NA[missing_idx_list[[j]]] <- NA

  # fit
  set.seed(10)
  init <- eSVD::initialization(dat_impute_NA, family = "curved_gaussian",
                               k = k, max_val = max_val)
  fit <- eSVD::fit_factorization(dat_impute_NA, u_mat = init$u_mat, v_mat = init$v_mat,
                                           family = "curved_gaussian",
                                           max_iter = 50, max_val = max_val,
                                           scalar = scalar_val,
                                           return_path = F, cores = ncores,
                                           verbose = T)
  pred_mat <- fit$u_mat %*% t(fit$v_mat)
  pred_val <- 1/pred_mat[missing_idx_list[[j]]]
  observed_val <- dat_impute[missing_idx_list[[j]]]

  res_list[[j]] <- cbind(pred_val, observed_val)

  save.image(paste0("../results/step3_scalar_heuristic", suffix, "_tmp.RData"))
}

rm(list = c("j", "init", "fit", "dat_impute_NA", "scalar_val_init"))
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




