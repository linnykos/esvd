set.seed(10)
load(paste0("../results/step3_scalar_heuristic", suffix, ".RData"))

init <- eSVD::initialization(dat_impute, family = family, k = k, max_val = max_val)
res_our <- eSVD::fit_factorization(dat_impute, u_mat = init$u_mat, v_mat = init$v_mat,
                                                 family = family,  reparameterize = T,
                                                 max_iter = 100, max_val = max_val,
                                                 scalar = scalar_val,
                                                 return_path = F, cores = ncores,
                                                 verbose = T)

rm(list = c("dat_impute_NA"))
save.image(paste0("../results/step4_factorization", suffix, ".RData"))
print(paste0(Sys.time(), ": Finished factorizing"))
