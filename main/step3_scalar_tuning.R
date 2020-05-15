set.seed(10)
load(paste0("../results/step2_rescaling", suffix, ".RData"))
session_info <- sessionInfo(); date_of_run <- Sys.time()

paramMat_esvd <- as.matrix(expand.grid(c(0.5, 1, 2, 4), c(3, 5, 10, 20)))
colnames(paramMat_esvd) <- c("scalar", "k")
esvd_missing_list <- vector("list", nrow(paramMat_esvd))
fitting_distr <- "curved_gaussian"

for(i in 1:nrow(paramMat_esvd)){
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
                                   return_path = F, cores = ncores,
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
