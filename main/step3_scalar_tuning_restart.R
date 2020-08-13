set.seed(10)
load(paste0("../results/step3_scalar_tuning", suffix, "_tmp.RData"))
ncores <- 5
start_idx <- max(which(sapply(esvd_missing_list, length) > 0)) + 1

for(i in start_idx:nrow(paramMat_esvd)){
  print(paste0("On parameter setting row ", i))
  tmp_list <- vector("list", length(cv_trials))

  for(j in 1:cv_trials){
    print(paste0("On trial ", j))

    # set missing values
    dat_impute_NA <- dat_impute
    dat_impute_NA[missing_idx_list[[j]]] <- NA

    # fit
    set.seed(10)
    init <- eSVD::initialization(dat_impute_NA, family = fitting_distr,
                                 k = paramMat_esvd[i,"k"], max_val = max_val, scalar = paramMat_esvd[i,"scalar"])
    tmp_list[[j]] <- eSVD::fit_factorization(dat_impute_NA, u_mat = init$u_mat, v_mat = init$v_mat,
                                             family = fitting_distr,
                                             max_iter = 50, max_val = max_val,
                                             scalar = paramMat_esvd[i,"scalar"],
                                             return_path = F, ncores = ncores,
                                             verbose = T)
    save.image(paste0("../results/step3_scalar_tuning", suffix, "_tmp.RData"))
  }

  esvd_missing_list[[i]] <- tmp_list
}


rm(list = c("j", "i", "init", "tmp_list", "dat_impute_NA"))
print(paste0(Sys.time(), ": Finished scalar heuristic"))
source_code_info <- c(source_code_info, readLines("../main/step3_scalar_tuning.R"))
save.image(paste0("../results/step3_scalar_tuning", suffix, ".RData"))
print(warnings())
