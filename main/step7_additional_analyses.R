# load(paste0("../results/step5_trajectory", suffix, ".RData"))
load(paste0("../results/step4_factorization", suffix, ".RData"))

dat_se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = t(dat_count)))
zinbwave_res <- zinbwave::zinbwave(dat_se, K = 5, maxiter.optimize = 100, normalizedValues = T,
                          commondispersion = F)
save.image(paste0("../results/step7_additional_analyses", suffix, ".RData"))

zinbwave_embedding <- SingleCellExperiment::reducedDims(zinbwave_res)$zinbwave
save.image(paste0("../results/step7_additional_analyses", suffix, ".RData"))
print(paste0(Sys.time(), ": Finished ZINB-WaVE"))

# include UMAPs here

print(warnings())
