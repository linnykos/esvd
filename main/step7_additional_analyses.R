load(paste0("../results/step5_trajectory", suffix, ".RData"))

dat_se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = t(dat_org)))
tmp <- zinbwave::zinbwave(dat_se, K = 5, maxiter.optimize = 100, normalizedValues = T,
                          commondispersion = F)
save.image("../experiment/Writeup_revision7_zinbwave.RData")
