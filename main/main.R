rm(list=ls())
set.seed(10)
library(eSVD)

suffix <- ""
ncores <- 20
doMC::registerDoMC(cores = ncores)

source("../main/step0_screening.R")
source("../main/step1_imputing.R")
source("../main/step2_naive_svd.R")
source("../main/step3_scalar_heuristic.R")
source("../main/step4_factorization.R")
source("../main/step5_clustering.R")
# source("../main/step6_figures.R")
