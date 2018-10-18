rm(list=ls())
source("../experiment/Week27_simulation_generator.R")
library(singlecell)

set.seed(10)
simulation <- .data_generator_exponential(total = 100, col_drop = F)
dat <- simulation$dat

zz <- dat[dat > 0]
max_val <- -1/mean(zz[zz < quantile(zz, probs = 0.2)])
max_iter <- 100

set.seed(10)
res_ideal <- singlecell:::.fit_factorization(simulation$obs_mat, simulation$cell_mat, simulation$gene_mat,
                                             verbose = T, family = "exponential",
                                             max_iter = max_iter, tol = NA,
                                             cores = NA, max_val = -max(abs(simulation$gram_mat)))

save.image("Week27_simulation.RData")

##########

set.seed(10)
init <- singlecell:::.initialization(dat, family = "exponential",
                                     max_val = max_val)
res_nodropout <- singlecell:::.fit_factorization(dat, init$u_mat, init$v_mat,
                                                 verbose = T, family = "exponential",
                                                 max_iter = max_iter, tol = NA,
                                                 cores = 15, max_val = max_val)

save.image("Week27_simulation.RData")

######################

set.seed(10)
dropout_mat <- singlecell:::.dropout(dat)
zero_mat <- singlecell:::.find_true_zeros(dropout_mat, num_neighbors = 14)

set.seed(10)
dat2 <- dat
dat2[which(is.na(zero_mat))] <- NA
dat_impute <- singlecell:::.scImpute(dat, which(is.na(zero_mat)), Kcluster = 4,
                                      verbose = T, weight = 0.25)
init_impute <- singlecell:::.initialization(dat_impute, family = "exponential",
                                            max_val = max_val)
save.image("Week27_simulation.RData")

set.seed(10)
res_withimpute <- singlecell:::.fit_factorization(dat_impute, init_impute$u_mat, init_impute$v_mat,
                                                  verbose = T, family = "exponential",
                                                  max_iter = max_iter, tol = NA,
                                                  cores = 15, max_val = max_val)

save.image("Week27_simulation.RData")

###################

dat_impute2 <- singlecell:::.scImpute(dat, which(is.na(zero_mat)), Kcluster = 4,
                                     verbose = T, weight = 0.5)
init_impute2 <- singlecell:::.initialization(dat_impute2, family = "exponential",
                                            max_val = max_val)
save.image("Week27_simulation.RData")

set.seed(10)
res_withimpute2 <- singlecell:::.fit_factorization(dat_impute2, init_impute2$u_mat, init_impute2$v_mat,
                                                  verbose = T, family = "exponential",
                                                  max_iter = max_iter, tol = NA,
                                                  cores = 15, max_val = max_val)

save.image("Week27_simulation.RData")
