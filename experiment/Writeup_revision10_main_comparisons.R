rm(list=ls())
load("../results/step4_factorization_original.RData")

zz = dat_impute
load("../results/step5_trajectory.RData")

dim(zz)
dim(dat_impute)

quantile(zz)
quantile(dat_impute)
