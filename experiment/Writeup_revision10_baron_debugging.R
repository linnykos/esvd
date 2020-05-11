rm(list=ls())
load("../results/step2_baron_scalar_tuning_spca_descend.RData")

ncores <- 10
doMC::registerDoMC(cores = ncores)

fitting_func <- function(dat_impute, vec, missing_idx_list){
  fitting_distr <- c("gaussian", "neg_binom", "curved_gaussian")[vec["fitting_distr"]]

  dat_org <- dat_impute*1000/max(dat_impute)
  dat_NA <- dat_org

  tmp_list <- vector("list", length(cv_trials))
  for(j in 2){
    dat_NA[missing_idx_list[[j]]] <- NA

    set.seed(10)
    init <- eSVD::initialization(dat_NA, family = fitting_distr, k = vec["k"], max_val = vec["max_val"],
                                 scalar = vec["scalar"])
    tmp_list[[j]] <- eSVD::fit_factorization(dat_NA, u_mat = init$u_mat, v_mat = init$v_mat,
                                             family = fitting_distr, max_iter = vec["max_iter"],
                                             scalar = vec["scalar"],
                                             max_val = vec["max_val"], return_path = F, cores = ncores, verbose = T)
  }

  tmp_list
}

i <- 3
ii <- 6
tmp[[ii]] <- fitting_func(dat_impute = preprocessing_list[[i]]$dat_impute, vec = paramMat_esvd[ii,],
                          missing_idx_list = missing_idx_list_list[[i]])

#############

i <- 3
ii <- 6
j <-3
dat_impute = preprocessing_list[[i]]$dat_impute
vec = paramMat_esvd[ii,]
missing_idx_list = missing_idx_list_list[[i]]

fitting_distr <- c("gaussian", "neg_binom", "curved_gaussian")[vec["fitting_distr"]]
dat_org <- dat_impute*1000/max(dat_impute)
dat_NA <- dat_org
dat_NA[missing_idx_list[[j]]] <- NA
set.seed(10)
init <- eSVD::initialization(dat_NA, family = fitting_distr, k = vec["k"], max_val = vec["max_val"],
                             scalar = vec["scalar"])
tmp <- eSVD::fit_factorization(dat_NA, u_mat = init$u_mat, v_mat = init$v_mat,
                                         family = fitting_distr, max_iter = vec["max_iter"],
                                         scalar = vec["scalar"],
                                         max_val = vec["max_val"], return_path = F, cores = ncores, verbose = T)
