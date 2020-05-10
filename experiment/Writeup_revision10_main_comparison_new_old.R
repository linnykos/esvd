rm(list=ls())
load("../results/step5_trajectory.RData")

ncores <- 20
doMC::registerDoMC(cores = ncores)

set.seed(10)
init <- eSVD::initialization(dat_impute, family = "curved_gaussian", k = 5, max_val = 5000,
                             scalar = 2)
esvd_embedding <- eSVD::fit_factorization(dat_impute, u_mat = init$u_mat, v_mat = init$v_mat,
                                          family = "curved_gaussian",
                                          max_iter = 100, max_val = 5000,
                                          scalar = 2,
                                          return_path = F, cores = ncores,
                                          verbose = T)

save.image("../results/Writeup_revision10_main_comparison_new_old.RData")
