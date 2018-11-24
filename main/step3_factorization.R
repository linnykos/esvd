load("../results/step1_imputing.RData")

dat_impute_log <- log(dat_impute + 1)

init <- singlecell:::.initialization(dat_impute_log, family = "gaussian", max_val = 10,
                                     k = 5)
res <- singlecell:::.fit_factorization(dat_impute_log, init$u_mat, init$v_mat,
                                       max_val = 10, family = "gaussian", verbose = T,
                                       max_iter = 25, reparameterize = T,
                                       return_path = F, cores = 10)

save.image("../results/step3_factorization_logged.RData")
