rm(list=ls())
source("../experiment/Week24_simulation_generator.R")
library(singlecell)

set.seed(10)
simulation <- .data_generator(total = 200, col_drop = F)
dat <- simulation$dat

init <- singlecell:::.initialization(dat, family = "gaussian")
res_nodropout <- singlecell:::.fit_factorization(dat, init$u_mat, init$v_mat,
                                                 verbose = T, family = "gaussian",
                                                 cores = 15)

save.image("Week24_simulation.RData")

dropout_mat <- singlecell:::.dropout(dat)
zero_mat <- singlecell:::.find_true_zeros(dropout_mat)

dat2 <- dat
dat2[which(is.na(zero_mat))] <- NA
dat_impute <- singlecell:::.scImpute(dat, which(is.na(zero_mat)), Kcluster = 4,
                                     max_time = 5, verbose = T)

save.image("Week24_simulation.RData")

init_impute <- singlecell:::.initialization(dat_impute, family = "gaussian")
res_withdropout <- singlecell:::.fit_factorization(dat2, init_impute$u_mat, init_impute$v_mat,
                                                   verbose = T, family = "gaussian",
                                                   cores = 15)

save.image("Week24_simulation.RData")

res_withimpute <- singlecell:::.fit_factorization(dat_impute, init_impute$u_mat, init_impute$v_mat,
                                                  verbose = T, family = "gaussian",
                                                  cores = 15)

save.image("Week24_simulation.RData")
