rm(list=ls())
set.seed(10)
library(eSVD)

suffix <- ""
ncores <- 5
doMC::registerDoMC(cores = ncores)

session_info <- sessionInfo()
source_code_info <- ""
date_of_run <- Sys.time()

set.seed(10)
load(paste0("../results/step3_scalar_tuning", suffix, "_tmp.RData"))
session_info <- sessionInfo(); date_of_run <- Sys.time()

scalar <- 2
k <- 5

set.seed(10)
init <- eSVD::initialization(dat_impute, family = fitting_distr, k = k, max_val = max_val,
                             scalar = scalar)
esvd_embedding <- eSVD::fit_factorization(dat_impute, u_mat = init$u_mat, v_mat = init$v_mat,
                                          family = fitting_distr,
                                          max_iter = 100, max_val = max_val,
                                          scalar = scalar,
                                          return_path = F, ncores = ncores,
                                          verbose = T)

source_code_info <- c(source_code_info, readLines("../main/step4_factorization.R"))
print(paste0(Sys.time(), ": Finished factorizing"))
save.image(paste0("../experiment/Writeup_revision15_factorization.RData"))
print(warnings())