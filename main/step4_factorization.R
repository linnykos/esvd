set.seed(10)
load("../results/step3_scalar_heuristic.RData")

init <- singlecell::initialization(dat_impute, family = family, extra_weight = extra_weight,
                                     k = k, max_val = max_val)
res_our <- singlecell::fit_factorization(dat_impute, u_mat = init$u_mat, v_mat = init$v_mat,
                                                 family = family,  reparameterize = T,
                                                 max_iter = 100, max_val = max_val,
                                                 scalar = scalar_val, extra_weight = extra_weight,
                                                 return_path = F, cores = 15,
                                                 verbose = T)

rm(list = c("dat_impute_NA"))
print(paste0(Sys.time(), ": Finished factorizing"))
save.image(paste0("../results/step4_factorization", suffix, ".RData"))
