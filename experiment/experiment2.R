rm(list=ls())

set.seed(10)
library(eSVD)

suffix <- ""
ncores <- 20
doMC::registerDoMC(cores = ncores)

load(paste0("../results/step5_trajectory", suffix, ".RData"))

esvd_sd_val <- eSVD::compute_curve_sd(esvd_curves, esvd_bootstrap_list, cores = ncores, verbose = T)
