rm(list=ls())
set.seed(10)
library(eSVD)
library(Seurat)

suffix <- ""
ncores <- 20
doMC::registerDoMC(cores = ncores)

session_info <- sessionInfo(); date_of_run <- Sys.time()
source_code_info <- ""

# source("../main_zeisel/step0_zeisel_preprocessing.R")
# source("../main_zeisel/step1_zeisel_gaussian_fitting.R")
source("../main_zeisel/step2_zeisel_scalar_tuning.R")
