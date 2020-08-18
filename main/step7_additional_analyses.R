load(paste0("../results/step6_cascade", suffix, ".RData"))
session_info <- sessionInfo(); date_of_run <- Sys.time()

dat_se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = t(dat_count)))
zinbwave_res <- zinbwave::zinbwave(dat_se, K = 5, maxiter.optimize = 100, normalizedValues = T,
                          commondispersion = F)
save.image(paste0("../results/step7_additional_analyses", suffix, ".RData"))

zinbwave_embedding <- SingleCellExperiment::reducedDims(zinbwave_res)$zinbwave
save.image(paste0("../results/step7_additional_analyses", suffix, ".RData"))
print(paste0(Sys.time(), ": Finished ZINB-WaVE"))

# include UMAPs here
set.seed(10)
config <- umap::umap.defaults
config$n_neighbors <- 30
config$verbose <- T
umap_all <- umap::umap(dat_impute, config = config)

##############

# include negative binomail fits
fitting_distr2 <- "neg_binom"
paramMat_esvd2 <- as.matrix(expand.grid(c(500, 2500, 10000), c(3, 5, 10)))
colnames(paramMat_esvd2) <- c("scalar", "k")
esvd_missing_list2 <- vector("list", nrow(paramMat_esvd2))

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
                                 k = paramMat_esvd2[i,"k"], max_val = max_val, scalar = paramMat_esvd2[i,"scalar"])
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

nat_mat_list_list <- lapply(1:nrow(paramMat_esvd2), function(i){
  lapply(1:cv_trials, function(j){
    u_mat <- esvd_missing_list2[[i]][[j]]$u_mat
    v_mat <- esvd_missing_list2[[i]][[j]]$v_mat
    u_mat %*% t(v_mat)
  })
})

esvd_angle_res2 <- eSVD::tuning_select_scalar(dat = dat_impute, nat_mat_list_list = nat_mat_list_list,
                                             family = fitting_distr2,  missing_idx_list = missing_idx_list,
                                             scalar_vec = paramMat_esvd2[,"scalar"])

rm(list = c("j", "i", "init", "tmp_list", "dat_impute_NA", "nat_mat_list_list"))
print(paste0(Sys.time(), ": Alternative analyses"))
source_code_info <- c(source_code_info, readLines("../main/step7_additional_analyses.R"))
save.image(paste0("../results/step7_additional_analyses", suffix, ".RData"))
print(warnings())

