rm(list=ls())
set.seed(10)
library(eSVD)

sessionInfo()

suffix <- ""
ncores <- 20
doMC::registerDoMC(cores = ncores)

source("../main_supplement/step0_baron_preprocessing.R")
