set.seed(10)
load(paste0("../results/step3_scalar_heuristic", suffix, ".RData"))

esvd_angle_vec <- lapply(1:nrow(paramMat_esvd), function(i){
  nat_mat_list <- lapply(1:cv_trials, function(j){
    u_mat <- esvd_missing_list[[i]][[j]]$u_mat
    v_mat <- esvd_missing_list[[i]][[j]]$v_mat
    u_mat %*% t(v_mat)
  })

  eSVD::plot_prediction_against_observed(dat_impute, nat_mat_list, family = fitting_distr,
                                         missing_idx_list, scalar = paramMat_esvd[i, "scalar"], plot = F)
})

esvd_angle_vec <- cbind(sapply(esvd_angle_vec, function(x){x$angle_val}),
                        sapply(esvd_angle_vec, function(x){x$bool}))

# REVISIT: for now, ignore the possibility of bool being false
idx <- which.min(abs(esvd_angle_vec[,1] - 45))
k <- paramMat_esvd[idx, "k"]
scalar <- paramMat_esvd[idx, "scalar"]

init <- eSVD::initialization(dat_impute, family = fitting_distr, k = k, max_val = max_val)
esvd_embedding <- eSVD::fit_factorization(dat_impute, u_mat = init$u_mat, v_mat = init$v_mat,
                                                 family = fitting_distr,
                                                 max_iter = 100, max_val = max_val,
                                                 scalar = scalar,
                                                 return_path = F, cores = ncores,
                                                 verbose = T)

rm(list = c("nat_mat_list", "idx", "init"))
save.image(paste0("../results/step4_factorization", suffix, ".RData"))
print(paste0(Sys.time(), ": Finished factorizing"))
