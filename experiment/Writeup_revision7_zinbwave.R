rm(list=ls())
load("../results/step5_trajectory.RData")
dat_se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = t(dat_impute)))
tmp <- zinbwave::zinbwave(dat_se, K = 5, maxiter.optimize = 100, normalizedValues = F,
                          commondispersion = F)
save("../experiment/Writeup_revision7_zinbwave.RData")
