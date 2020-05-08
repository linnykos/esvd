# load(paste0("../results/step5_trajectory", suffix, ".RData"))
load(paste0("../results/step4_factorization", suffix, ".RData"))

dat_se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = t(dat_count)))
zinbwave_res <- zinbwave::zinbwave(dat_se, K = 5, maxiter.optimize = 100, normalizedValues = T,
                          commondispersion = F)
save.image("../experiment/Writeup_revision7_zinbwave.RData")

zinbwave_embedding <- SingleCellExperiment::reducedDims(zinbwave_res)$zinbwave
save.image("../experiment/Writeup_revision7_zinbwave.RData")

# include UMAPs here
