set.seed(10)
load(paste0("../results/step4_factorization_original.RData"))

suffix <- ""
ncores <- 20
doMC::registerDoMC(cores = ncores)
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

rm(list = c("nat_mat_list_list", "idx"))
source_code_info <- c(source_code_info, readLines("../main/step4_factorization.R"))
print(paste0(Sys.time(), ": Finished factorizing"))
save.image(paste0("../experiment/Writeup_revision15_debugging_version.RData"))
print(warnings())
