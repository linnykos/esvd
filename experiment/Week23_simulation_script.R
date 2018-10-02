rm(list=ls())
source("../experiment/Week23_simulation_generator.R")

library(singlecell)
set.seed(10)
simulation <- .data_generator(total = 200)
dat <- simulation$dat

init <- singlecell:::.initialization(dat)
res_nodropout <- singlecell:::.fit_exponential_factorization(dat, init$u_mat, init$v_mat,
                                                             verbose = T)

save.image("Week23_simulation_2.RData")

dropout_mat <- singlecell:::.dropout(dat)
zero_mat <- singlecell:::.find_true_zeros(dropout_mat)

dat2 <- dat
dat2[which(is.na(zero_mat))] <- NA
dat_impute <- singlecell:::.scImpute(dat, which(is.na(zero_mat)), Kcluster = 4,
                                     max_time = 5)

init_impute <- singlecell:::.initialization(dat_impute)
res_withdropout <- singlecell:::.fit_exponential_factorization(dat2, init_impute$u_mat,
                                                               init_impute$v_mat,
                                                               verbose = T)

save.image("Week23_simulation_2.RData")
