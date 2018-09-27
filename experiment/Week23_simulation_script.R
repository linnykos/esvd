rm(list=ls())
source("../experiment/Week23_simulation_generator.R")
library(singlecell)

set.seed(10)
simulation <- .data_generator(total = 150, distr_func = function(x){stats::rnorm(1, x, x/2)})
dat <- simulation$dat

res_nodropout <- singlecell:::.fit_gaussian_factorization(dat, k = 5, verbose = T)

dropout_mat <- singlecell:::.dropout(dat)
zero_mat <- singlecell:::.find_true_zeros(dropout_mat)

dat2 <- dat
dat2[which(is.na(zero_mat))] <- NA
res_withdropout <- singlecell:::.fit_gaussian_factorization(dat2, k = 5, verbose = T)

save.image("../experiment/Week23_tmp.RData")
