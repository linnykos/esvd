rm(list=ls())
set.seed(10)
library(eSVD)
library(Seurat)

sessionInfo()

suffix <- "_original"
ncores <- 20
doMC::registerDoMC(cores = ncores)

# source("../main/step0_screening.R")
# source("../main/step1_naive_svd.R")
# source("../main/step2_rescaling.R")
# source("../main/step3_scalar_tuning.R")
# source("../main/step4_factorization.R")
# source("../main/step5_trajectory.R")
# source("../main/step6_figures.R")
source("../main/step7_additional_analyses.R")
