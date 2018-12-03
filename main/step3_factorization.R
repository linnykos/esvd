set.seed(10)
load("../results/step2b_scalar_heuristic.RData")

init <- singlecell:::.initialization(dat_impute, family = "gaussian", scalar = scalar_vec[i],
                                     k = k, max_val = max_val)
res_our <- singlecell:::.fit_factorization(dat_impute, u_mat = init$u_mat, v_mat = init$v_mat,
                                                 family = "gaussian",  reparameterize = T,
                                                 max_iter = 100, max_val = max_val,
                                                 scalar = scalar_val,
                                                 return_path = F, cores = 15,
                                                 verbose = T)

save.image("../results/step3_factorization.RData")
