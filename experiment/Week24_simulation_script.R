rm(list=ls())
source("../experiment/Week24_simulation_generator.R")
library(singlecell)

set.seed(10)
simulation <- .data_generator_exponential(total = 200, col_drop = F)
dat <- simulation$dat

init_ideal <- singlecell:::.initialization(simulation$obs_mat, family = "exponential")
res_ideal <- singlecell:::.fit_factorization(simulation$obs_mat, init_ideal$u_mat, init_ideal$v_mat,
                                             verbose = T, family = "exponential",
                                             cores = 15)

save.image("Week24_simulation_exponential.RData")

init <- singlecell:::.initialization(dat, family = "exponential")
res_nodropout <- singlecell:::.fit_factorization(dat, init$u_mat, init$v_mat,
                                                 verbose = T, family = "exponential",
                                                 cores = 15)

save.image("Week24_simulation_exponential.RData")

dropout_mat <- singlecell:::.dropout(dat)
zero_mat <- singlecell:::.find_true_zeros(dropout_mat, num_neighbors = 15)

dat2 <- dat
dat2[which(is.na(zero_mat))] <- NA
dat_impute <- singlecell:::.scImpute(dat, which(is.na(zero_mat)), Kcluster = 4,
                                     max_time = 5, verbose = T)

save.image("Week24_simulation_exponential.RData")

init_impute <- singlecell:::.initialization(dat_impute, family = "exponential")
res_withdropout <- singlecell:::.fit_factorization(dat2, init_impute$u_mat, init_impute$v_mat,
                                                   verbose = T, family = "exponential",
                                                   cores = 15)

save.image("Week24_simulation_exponential.RData")

res_withimpute <- singlecell:::.fit_factorization(dat_impute, init_impute$u_mat, init_impute$v_mat,
                                                  verbose = T, family = "exponential",
                                                  cores = 15)

save.image("Week24_simulation_exponential.RData")

