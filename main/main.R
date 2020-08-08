rm(list=ls())
set.seed(10)
library(eSVD)

suffix <- "_debugging2"
ncores <- 20
doMC::registerDoMC(cores = ncores)

session_info <- sessionInfo()
source_code_info <- ""
date_of_run <- Sys.time()

source("../main/step4_factorization.R")
# source("../main/step5_trajectory.R")
# source("../main/step6_additional_analyses.R")
