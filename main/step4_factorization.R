set.seed(10)
load(paste0("../results/step3_scalar_tuning", suffix, ".RData"))

nat_mat_list_list <- lapply(1:nrow(paramMat_esvd), function(i){
  lapply(1:cv_trials, function(j){
    u_mat <- esvd_missing_list[[i]][[j]]$u_mat
    v_mat <- esvd_missing_list[[i]][[j]]$v_mat
    u_mat %*% t(v_mat)
  })
})

esvd_angle_res <- eSVD::tuning_select_scalar(dat = dat_impute, nat_mat_list_list = nat_mat_list_list,
                                family = fitting_distr,  missing_idx_list = missing_idx_list,
                                scalar_vec = paramMat_esvd[,"scalar"])
scalar <- paramMat_esvd[esvd_angle_res$idx, "scalar"]
k <- paramMat_esvd[esvd_angle_res$idx, "k"]

init <- eSVD::initialization(dat_impute, family = fitting_distr, k = k, max_val = max_val,
                             scalar = scalar)
esvd_embedding <- eSVD::fit_factorization(dat_impute, u_mat = init$u_mat, v_mat = init$v_mat,
                                                 family = fitting_distr,
                                                 max_iter = 100, max_val = max_val,
                                                 scalar = scalar,
                                                 return_path = F, cores = ncores,
                                                 verbose = T, tol = 1e-4)

rm(list = c("nat_mat_list_list", "idx"))
save.image(paste0("../results/step4_factorization", suffix, ".RData"))
print(paste0(Sys.time(), ": Finished factorizing"))
print(warnings())

##########################
#
# i <- 7
# nat_mat_list <- lapply(1:cv_trials, function(j){
#   u_mat <- esvd_missing_list[[i]][[j]]$u_mat
#   v_mat <- esvd_missing_list[[i]][[j]]$v_mat
#   u_mat %*% t(v_mat)
# })
# plot_prediction_against_observed(dat_impute, nat_mat_list = nat_mat_list,
#                                  family = fitting_distr, missing_idx_list = missing_idx_list,
#                                  scalar = paramMat_esvd[i, "scalar"])

