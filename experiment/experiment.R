rm(list=ls())
load("../results/step2_naive_svd_spca.RData")

max_val <- 5000
scalar_vec <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 100)
res_list <- vector("list", length(scalar_vec))

i <- 1
print(paste0("Trying scalar value = ", scalar_vec[i]))
init <- singlecell::initialization(dat_impute_NA, family = family,
                                   k = k, max_val = max_val)

stopifnot(all(init$u_mat %*% t(init$v_mat) > 0))
range(init$u_mat %*% t(init$v_mat) )

res_list[[i]] <- singlecell::fit_factorization(dat_impute_NA, u_mat = init$u_mat, v_mat = init$v_mat,
                                               family = family, reparameterize = T,
                                               max_iter = 25, max_val = NA,
                                               scalar = scalar_vec[i],
                                               return_path = F, cores = NA,
                                               verbose = T)
##########################

