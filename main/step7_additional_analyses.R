load(paste0("../results/step5_trajectory", suffix, ".RData"))

dat_se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = t(dat_count)))
zinbwave_res <- zinbwave::zinbwave(dat_se, K = 5, maxiter.optimize = 100, normalizedValues = T,
                          commondispersion = F)
save.image(paste0("../results/step7_additional_analyses", suffix, ".RData"))

zinbwave_embedding <- SingleCellExperiment::reducedDims(zinbwave_res)$zinbwave
save.image(paste0("../results/step7_additional_analyses", suffix, ".RData"))
print(paste0(Sys.time(), ": Finished ZINB-WaVE"))

# include UMAPs here

# include negative binomail fits
fitting_distr2 <- "neg_binom"
paramMat_esvd2 <- as.matrix(expand.grid(c(500, 2500, 10000), c(3, 5, 10)))
esvd_missing_list2 <- vector("list", nrow(paramMat_esvd))

for(i in 1:nrow(paramMat_esvd2)){
  print(paste0("On parameter setting row ", i))
  tmp_list <- vector("list", length(cv_trials))

  for(j in 1:cv_trials){
    print(paste0("On trial ", j))

    # set missing values
    dat_impute_NA <- dat_impute
    dat_impute_NA[missing_idx_list[[j]]] <- NA

    # fit
    set.seed(10)
    init <- eSVD::initialization(dat_impute_NA, family = fitting_distr2,
                                 k = paramMat_esvd[i,"k"], max_val = max_val, scalar = paramMat_esvd2[i,"scalar"])
    tmp_list[[j]] <- eSVD::fit_factorization(dat_impute_NA, u_mat = init$u_mat, v_mat = init$v_mat,
                                             family = fitting_distr2,
                                             max_iter = 50, max_val = max_val,
                                             scalar = paramMat_esvd2[i,"scalar"],
                                             return_path = F, cores = ncores,
                                             verbose = T)
    save.image(paste0("../results/step7_additional_analyses", suffix, "_tmp.RData"))
  }

  esvd_missing_list2[[i]] <- tmp_list
  save.image(paste0("../results/step7_additional_analyses", suffix, "_tmp.RData"))
}

rm(list = c("j", "i", "init", "tmp_list", "dat_impute_NA"))
print(paste0(Sys.time(), ": Alternative analyses"))
source_code_info <- c(source_code_info, readLines("../main/step7_additional_analyses.R"))
save.image(paste0("../results/step7_additional_analyses", suffix, ".RData"))
print(warnings())

