set.seed(10)
load(paste0("../results/step3_scalar_heuristic", suffix, ".RData"))

esvd_angle_vec <- rep(NA, nrow(paramMat))
for(i in 1:nrow(paramMat)){
  nat_mat_list <- lapply(1:cv_trials, function(j){
    u_mat <- esvd_missing_list[[i]][[j]]$u_mat
    v_mat <- esvd_missing_list[[i]][[j]]$v_mat
    u_mat %*% t(v_mat)
  })

  esvd_angle_vec[i] <- eSVD::plot_prediction_against_observed(dat_impute, nat_mat_list, family = "curved_gaussian",
                                         missing_idx_list, scalar = paramMat[i, "scalar"], plot = F)
}

idx <- which.min(abs(esvd_angle_vec - 45))
k <- paramMat[idx, "k"]
scalar <- paramMat[idx, "scalar"]

init <- eSVD::initialization(dat_impute, family = "curved_gaussian", k = k, max_val = max_val)
esvd_embedding <- eSVD::fit_factorization(dat_impute, u_mat = init$u_mat, v_mat = init$v_mat,
                                                 family = "curved_gaussian",
                                                 max_iter = 100, max_val = max_val,
                                                 scalar = scalar,
                                                 return_path = F, cores = ncores,
                                                 verbose = T)

rm(list = c("nat_mat_list", "idx", "init"))
save.image(paste0("../results/step4_factorization", suffix, ".RData"))
print(paste0(Sys.time(), ": Finished factorizing"))
