set.seed(10)
load(paste0("../results/step1_baron_gaussian_fitting", suffix, ".RData"))

neg_binom_vec <- c(250, 1000, 1e4, 1e6)
curved_gaussian_vec <- c(1, 2, 4, 100)

paramMat <- rbind(as.matrix(expand.grid(2, k_vec, neg_binom_vec, 50, 3000)),
                  as.matrix(expand.grid(3, k_vec, curved_gaussian_vec, 50, 3000)))
colnames(paramMat) <- c("fitting_distr", "k", "scalar_val", "max_iter", "max_val")

esvd_missing_list_list <- vector("list", length(preprocessing_list))

######################################

fitting_func <- function(dat_impute, vec, missing_idx_list){
  fitting_distr <- c("gaussian", "neg_binom", "curved_gaussian")[vec["fitting_distr"]]

  dat_org <- dat_impute*1000/max(dat_impute)
  dat_NA <- dat_org

  tmp_list <- vector("list", length(cv_trials))
  for(j in 1:cv_trials){
    dat_NA[missing_idx_list[[j]]] <- NA

    set.seed(10)
    init <- eSVD::initialization(dat_NA, family = fitting_distr, k = vec["k"], max_val = vec["max_val"],
                                 scalar = vec["scalar_val"])
    tmp_list[[j]] <- eSVD::fit_factorization(dat_NA, u_mat = init$u_mat, v_mat = init$v_mat,
                                             family = fitting_distr, max_iter = vec["max_iter"],
                                             scalar = vec["scalar_val"],
                                             max_val = vec["max_val"], return_path = F, cores = ncores, verbose = F)
  }

  tmp_list
}

for(i in 1:length(preprocessing_list)){
  print(paste0(Sys.time(), ": Starting dataset ", i))
  tmp <- vector("list", nrow(paramMat))

  for(ii in 1:nrow(paramMat)){
    print(paste0("Starting parameter row ", ii))
    tmp[[ii]] <- fitting_func(dat_impute = preprocessing_list[[i]]$dat_impute, vec = paramMat[ii,],
                              missing_idx_list = missing_idx_list_list[[i]])
  }

  esvd_missing_list_list[[i]] <- tmp
  save.image(paste0("../results/step2_baron_scalar_tuning", suffix, ".RData"))
}

print(paste0(Sys.time(), ": Finished scalar tuning for eSVD"))
save.image(paste0("../results/step2_baron_scalar_tuning", suffix, ".RData"))
print(warnings())

