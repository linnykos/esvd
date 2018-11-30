load("../results/step1_imputing.RData")

#dat_impute_log <- log(dat_impute + 1)

max_val <- 15
init <- singlecell:::.initialization(dat_impute, family = "gaussian", max_val = max_val,
                                     k = 5)
res <- singlecell:::.fit_factorization(dat_impute, init$u_mat, init$v_mat,
                                       max_val = max_val, family = "gaussian", verbose = T,
                                       max_iter = 25, reparameterize = T,
                                       return_path = F, cores = 10)

save.image("../results/step3_factorization_logged.RData")
