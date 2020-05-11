rm(list=ls())
set.seed(10)
library(eSVD)
library(Seurat)

sessionInfo()

suffix <- "_spca_descend"
ncores <- 10
doMC::registerDoMC(cores = ncores)

session_info <- sessionInfo()
source_code_info <- ""
date_of_run <- Sys.time()

source("../main_supplement/step0_baron_preprocessing.R")
source("../main_supplement/step1_baron_gaussian_fitting.R")
source("../main_supplement/step2_baron_scalar_tuning.R")
