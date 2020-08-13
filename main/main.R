rm(list=ls())
set.seed(10)
library(eSVD)

suffix <- ""
ncores <- 5
doMC::registerDoMC(cores = ncores)

session_info <- sessionInfo()
source_code_info <- ""
date_of_run <- Sys.time()

# source("../main/step0_screening.R")
# source("../main/step1_naive_svd.R")
# source("../main/step2_rescaling.R")
# source("../main/step3_scalar_tuning.R")
source("../main/step3_scalar_tuning_restart.R")
source("../main/step4_factorization.R")
source("../main/step5_trajectory.R")
source("../main/step6_cascade.R")
