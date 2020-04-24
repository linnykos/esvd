rm(list=ls())
load("../results/step3_scalar_heuristic_cg_hvg_tmp.RData")

##############################

training_idx_list <- lapply(1:cv_trials, function(i){
  c(1:prod(dim(dat_impute)))[-missing_idx_list[[i]]]
})

nat_mat_list_list <- lapply(1:nrow(paramMat_esvd), function(i){
  lapply(1:cv_trials, function(j){
    u_mat <- esvd_missing_list[[i]][[j]]$u_mat
    v_mat <- esvd_missing_list[[i]][[j]]$v_mat
    u_mat %*% t(v_mat)
  })
})

esvd_angle_res <- eSVD:::tuning_select_scalar(dat = dat_impute, nat_mat_list_list = nat_mat_list_list,
                                              family = fitting_distr,  missing_idx_list = training_idx_list,
                                              scalar_vec = paramMat_esvd[,"scalar"])
# okay good, the training results are also pretty ass
