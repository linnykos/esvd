set.seed(10)
load(paste0("../results/step2_zeisel_scalar_tuning", suffix, ".RData"))
session_info <- sessionInfo(); date_of_run <- Sys.time()

dat <- dat * 1000/max(dat)

nat_mat_list_list <- lapply(1:nrow(paramMat_esvd), function(i){
  lapply(1:cv_trials, function(j){
    u_mat <- esvd_missing_list[[i]][[j]]$u_mat
    v_mat <- esvd_missing_list[[i]][[j]]$v_mat
    u_mat %*% t(v_mat)
  })
})

training_idx_list <- lapply(1:length(missing_idx_list), function(i){
  c(1:prod(dim(dat)))[-missing_idx_list[[i]]]
})

selection_mat <- as.matrix(expand.grid(k_vec, c(2,3)))
tuning_idx <- rep(NA, nrow(selection_mat))
for(j in 1:nrow(selection_mat)){
  distr_num <- selection_mat[j,2]
  k <- selection_mat[j,1]
  param_idx <- intersect(which(paramMat_esvd[,"fitting_distr"] == distr_num), which(paramMat_esvd[,"k"] == k))

  fitting_distr <- c("gaussian", "neg_binom", "curved_gaussian")[distr_num]
  esvd_angle_res <- eSVD::tuning_select_scalar(dat = dat, nat_mat_list_list = nat_mat_list_list[param_idx],
                                               family = fitting_distr,  missing_idx_list = missing_idx_list,
                                               scalar_vec = paramMat_esvd[param_idx,"scalar"])
  tuning_idx[j] <- esvd_angle_res$idx
}

# now plot
# par(mfrow = c(2,3))
# for(j in 1:nrow(selection_mat)){
#   print(j)
#
#   distr_num <- selection_mat[j,2]
#   k <- selection_mat[j,1]
#   param_idx <- intersect(which(paramMat_esvd[,"fitting_distr"] == distr_num), which(paramMat_esvd[,"k"] == k))
#   fitting_distr <- c("gaussian", "neg_binom", "curved_gaussian")[distr_num]
#   if(fitting_distr == "neg_binom"){
#     scalar <- neg_binom_vec[tuning_idx[j]]
#   } else {
#     scalar <- curved_gaussian_vec[tuning_idx[j]]
#   }
#   main_str <- paste0(ifelse(fitting_distr == "neg_binom", "Negative binomial", "Curved Gaussian"),
#                      "\n(k = ", k, ", scalar = ", scalar, ")")
#
#   eSVD::plot_prediction_against_observed(dat, nat_mat_list = nat_mat_list_list[param_idx][[tuning_idx[j]]],
#                                    missing_idx_list = missing_idx_list,
#                                    family = fitting_distr, scalar = scalar,
#                                    main = main_str,
#                                    max_points = 1e6)
# }

fitting_vec <- paramMat_esvd[12,]
fitting_distr <- c("gaussian", "neg_binom", "curved_gaussian")[fitting_vec["fitting_distr"]]
dat <- dat*1000/max(dat)

set.seed(10)
init <- eSVD::initialization(dat, family = fitting_distr, k = k,
                             max_val = fitting_vec["max_val"],
                             scalar = fitting_vec["scalar"])
esvd_embedding <- eSVD::fit_factorization(dat, u_mat = init$u_mat, v_mat = init$v_mat,
                                                    family = fitting_distr, max_iter = fitting_vec["max_iter"],
                                                    scalar = fitting_vec["scalar"],
                                                    max_val = fitting_vec["max_val"],
                                                    return_path = F, ncores = ncores, verbose = T)


rm("init", "fitting_distr", "nat_mat_list_list", "training_idx_list", "esvd_angle_res")
print(paste0(Sys.time(), ": Finished fitting eSVD"))
source_code_info <- c(source_code_info, readLines("../main_zeisel/step3_zeisel_factorization.R"))
save.image(paste0("../results/step3_zeisel_factorization", suffix, ".RData"))
print(warnings())


