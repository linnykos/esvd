load("../results/step1_imputing.RData")

init <- singlecell:::.initialization(dat_impute, family = "gaussian", max_val = 10,
                                     k = 5)
res <- singlecell:::.fit_factorization(dat_impute, init$u_mat, init$v_mat,
                                       max_val = 10, family = "gaussian", verbose = T,
                                       max_iter = 50, reparameterize = T,
                                       return_path = F)

save.image("../results/step3_factorization.RData")
