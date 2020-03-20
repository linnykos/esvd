rm(list=ls())

load("../results/lingxue_analysis.RData")

sapply(fit_all_list, length)
sapply(fit_list, length)
sapply(fit_list[[1]]$neg_binom_missing, length)
