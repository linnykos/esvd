rm(list=ls())
set.seed(10)
library(eSVD)

sessionInfo()

suffix <- "_cg_spca-vst_before_rescaling_300_all"
ncores <- 20
doMC::registerDoMC(cores = ncores)

source("../main/step0_screening.R")
source("../main/step1_naive_svd.R")
source("../main/step2_rescaling.R")
source("../main/step3_scalar_heuristic.R")
source("../main/step4_factorization.R")
# source("../main/step5_trajectory.R")
# source("../main/step6_figures.R")
