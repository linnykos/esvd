set.seed(10)
load(paste0("../results/step2_baron_scalar_tuning", suffix, ".RData"))
session_info <- sessionInfo(); date_of_run <- Sys.time()

# hand pick the winning factorizations
selected_paramMat <- cbind(c(3,  3, 3, 3, 3, 3),
                           c(10, 3, 3, 5, 5, 3),
                           c(2,  2, 4, 4, 4, 4), 50, 3000)
colnames(selected_paramMat) <- c("fitting_distr", "k", "scalar", "max_iter", "max_val")

esvd_embedding_list <- vector("list", length(preprocessing_list))
for(i in 1:length(preprocessing_list)){
  print(paste0(Sys.time(), ": Starting dataset ", i))

  fitting_distr <- c("gaussian", "neg_binom", "curved_gaussian")[vec["fitting_distr"]]
  dat_impute <- preprocessing_list[[i]]$dat_impute
  dat_impute <- dat_impute*1000/max(dat_impute)

  set.seed(10)
  init <- eSVD::initialization(dat_impute, family = fitting_distr, k = selected_paramMat[i, "k"],
                               max_val = selected_paramMat[i, "max_val"],
                               scalar = selected_paramMat[i, "scalar"])
  esvd_embedding_list[[i]] <- eSVD::fit_factorization(dat_impute, u_mat = init$u_mat, v_mat = init$v_mat,
                                           family = fitting_distr, max_iter = selected_paramMat[i, "max_iter"],
                                           scalar = selected_paramMat[i, "scalar"],
                                           max_val = selected_paramMat[i, "max_val"],
                                           return_path = F, cores = ncores, verbose = T)

  save.image(paste0("../results/step3_baron_factorization.RData"))
}


print(paste0(Sys.time(), ": Finished eSVD factorization"))
source_code_info <- c(source_code_info, readLines("../main_supplement/step3_baron_factorization.R"))
save.image(paste0("../results/step3_baron_factorization.RData"))
print(warnings())
