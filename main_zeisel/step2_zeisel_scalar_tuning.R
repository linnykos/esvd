set.seed(10)
load(paste0("../results/step1_zeisel_gaussian_fitting", suffix, ".RData"))
session_info <- sessionInfo(); date_of_run <- Sys.time()

neg_binom_vec <- c(250, 1000, 1e4, 1e6)
curved_gaussian_vec <- c(1, 2, 4, 100)

# paramMat_esvd <- rbind(as.matrix(expand.grid(2, k_vec, neg_binom_vec, 50, 3000)),
#                        as.matrix(expand.grid(3, k_vec, curved_gaussian_vec, 50, 3000)))
paramMat_esvd <- as.matrix(expand.grid(2, c(10, 20, 50), c(250, 500, 1000, 1e4, 5e4), 50, 3000))
colnames(paramMat_esvd) <- c("fitting_distr", "k", "scalar", "max_iter", "max_val")

######################################

fitting_func <- function(dat_impute, vec, missing_idx_list){
  fitting_distr <- c("gaussian", "neg_binom", "curved_gaussian")[vec["fitting_distr"]]

  dat_org <- dat_impute*1000/max(dat_impute)

  tmp_list <- vector("list", length(cv_trials))
  for(j in 1:cv_trials){

    # set missing value
    dat_NA <- dat_org
    print(paste0("Starting trial ", j))
    dat_NA[missing_idx_list[[j]]] <- NA

    set.seed(10)
    init <- eSVD::initialization(dat_NA, family = fitting_distr, k = vec["k"], max_val = vec["max_val"],
                                 scalar = vec["scalar"])
    tmp_list[[j]] <- eSVD::fit_factorization(dat_NA, u_mat = init$u_mat, v_mat = init$v_mat,
                                             family = fitting_distr, max_iter = vec["max_iter"],
                                             scalar = vec["scalar"],
                                             max_val = vec["max_val"], return_path = F, ncores = ncores, verbose = T)
  }

  tmp_list
}

esvd_missing_list <- vector("list", nrow(paramMat_esvd))

for(ii in 1:nrow(paramMat_esvd)){
  print(paste0("Starting parameter row ", ii))
  esvd_missing_list[[ii]] <- fitting_func(dat_impute = dat, vec = paramMat_esvd[ii,],
                            missing_idx_list = missing_idx_list)

  save.image(paste0("../results/step2_zeisel_scalar_tuning", suffix, ".RData"))
}

print(paste0(Sys.time(), ": Finished scalar tuning for eSVD"))
source_code_info <- c(source_code_info, readLines("../main_zeisel/step2_zeisel_scalar_tuning.R"))
save.image(paste0("../results/step2_zeisel_scalar_tuning", suffix, ".RData"))
print(warnings())

