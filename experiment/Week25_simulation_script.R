rm(list=ls())
source("../experiment/Week25_simulation_generator.R")
library(singlecell)
load("../experiment/Week25_simulation_exponential.RData")
#
# set.seed(10)
# simulation <- .data_generator_exponential(total = 200, col_drop = F)
# dat <- simulation$dat
#
# zz <- dat[dat > 0]
# max_val <- -1/mean(zz[zz < quantile(zz, probs = 0.2)])
# max_iter <- 50
#
# set.seed(10)
# res_ideal <- singlecell:::.fit_factorization(simulation$obs_mat, simulation$cell_mat, simulation$gene_mat,
#                                              verbose = T, family = "exponential",
#                                              max_iter = max_iter, tol = NA,
#                                              cores = NA, max_val = -max(abs(simulation$gram_mat)))
#
# save.image("Week25_simulation_exponential.RData")
#
# set.seed(10)
# init <- singlecell:::.initialization(dat, family = "exponential",
#                                      max_val = max_val)
# res_nodropout <- singlecell:::.fit_factorization(dat, init$u_mat, init$v_mat,
#                                                  verbose = T, family = "exponential",
#                                                  max_iter = max_iter, tol = NA,
#                                                  cores = 15, max_val = max_val)
#
# save.image("Week25_simulation_exponential.RData")
#
# res_nodropout_cheat <- singlecell:::.fit_factorization(dat, simulation$cell_mat, simulation$gene_mat,
#                                                  verbose = T, family = "exponential",
#                                                  max_iter = max_iter, tol = NA,
#                                                  cores = 15, max_val = max_val)
#
# save.image("Week25_simulation_exponential.RData")
#
# ######################
#
set.seed(10)
dropout_mat <- singlecell:::.dropout(dat)
zero_mat <- singlecell:::.find_true_zeros(dropout_mat, num_neighbors = 14)

set.seed(10)
dat2 <- dat
dat2[which(is.na(zero_mat))] <- NA
dat_impute <- singlecell:::.scImpute(dat, which(is.na(zero_mat)), Kcluster = 4,
                                     max_time = 5, verbose = T)

save.image("Week25_simulation_exponential.RData")

idx_impute2 <- which(is.na(zero_mat))
dat3 <- dat
dat[idx_impute2] <- -1/simulation$gram_mat
set.seed(10)
res_withimpute_population <- singlecell:::.fit_factorization(dat3, simulation$cell_mat, simulation$gene_mat,
                                                   verbose = T, family = "exponential",
                                                   max_iter = max_iter, tol = NA,
                                                   cores = 15, max_val = max_val)
save.image("Week25_simulation_exponential.RData")

# set.seed(10)
# init_impute <- singlecell:::.initialization(dat_impute, family = "exponential",
#                                             max_val = max_val)
# res_withdropout <- singlecell:::.fit_factorization(dat2, init_impute$u_mat, init_impute$v_mat,
#                                                    verbose = T, family = "exponential",
#                                                    max_iter = max_iter, tol = NA,
#                                                    cores = 15, max_val = max_val)
# save.image("Week25_simulation_exponential.RData")
#
# set.seed(10)
# res_withdropout_cheat <- singlecell:::.fit_factorization(dat2, simulation$cell_mat, simulation$gene_mat,
#                                                    verbose = T, family = "exponential",
#                                                    max_iter = max_iter, tol = NA,
#                                                    cores = 15, max_val = max_val)
#
# save.image("Week25_simulation_exponential.RData")
#
# ######
#
# set.seed(10)
# res_withimpute <- singlecell:::.fit_factorization(dat_impute, init_impute$u_mat, init_impute$v_mat,
#                                                   verbose = T, family = "exponential",
#                                                   max_iter = max_iter, tol = NA,
#                                                   cores = 15, max_val = max_val)
#
# save.image("Week25_simulation_exponential.RData")
#
# set.seed(10)
# res_withimpute_cheat <- singlecell:::.fit_factorization(dat_impute, simulation$cell_mat, simulation$gene_mat,
#                                                   verbose = T, family = "exponential",
#                                                   max_iter = max_iter, tol = NA,
#                                                   cores = 15, max_val = max_val)
#
# save.image("Week25_simulation_exponential.RData")
#
