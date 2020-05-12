rm(list=ls())
set.seed(10)
library(eSVD)

sessionInfo()

date_of_run <- Sys.time()
suffix <- "_original"
ncores <- 20
doMC::registerDoMC(cores = ncores)

session_info <- sessionInfo()
source_code_info <- ""
date_of_run <- Sys.time()

set.seed(10)
# load(paste0("../results/step3_scalar_tuning", suffix, ".RData"))
load("../results/step5_trajectory.RData")

nat_mat_list_list <- lapply(1:nrow(paramMat_esvd), function(i){
  lapply(1:cv_trials, function(j){
    u_mat <- esvd_missing_list[[i]][[j]]$u_mat
    v_mat <- esvd_missing_list[[i]][[j]]$v_mat
    u_mat %*% t(v_mat)
  })
})

fitting_distr <- "curved_gaussian"
esvd_angle_res <- eSVD::tuning_select_scalar(dat = dat_impute, nat_mat_list_list = nat_mat_list_list,
                                             family = fitting_distr,  missing_idx_list = missing_idx_list,
                                             scalar_vec = paramMat_esvd[,"scalar"])
scalar <- paramMat_esvd[esvd_angle_res$idx, "scalar"]
k <- paramMat_esvd[esvd_angle_res$idx, "k"]

set.seed(10)
init <- eSVD::initialization(dat_impute, family = fitting_distr, k = k, max_val = max_val,
                             scalar = scalar)
esvd_embedding <- eSVD::fit_factorization(dat_impute, u_mat = init$u_mat, v_mat = init$v_mat,
                                          family = fitting_distr,
                                          max_iter = 100, max_val = max_val,
                                          scalar = scalar,
                                          return_path = F, cores = ncores,
                                          verbose = T, tol = 1e-4)

rm(list = c("nat_mat_list_list", "idx"))
source_code_info <- c(source_code_info, readLines("../experiment/Writeup_revision10_main_comparison_new_old.R"))
save.image("../results/tmp.RData")
