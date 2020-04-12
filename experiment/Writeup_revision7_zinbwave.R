rm(list=ls())
load("../results/step0_screening.RData")

res_hvg <- descend::findHVG(res_descend, threshold = 50)
length(res_hvg$HVG.genes)

spca_summary <- cbind(v_seq, t(sapply(res_list, function(x){c(length(unique(sort(unlist(apply(x$v, 2, function(y){which(y != 0)}))))), x$prop.var.explained[5])})))
idx <- min(which(spca_summary[,2] == max(spca_summary[,2])))
target_var <- spca_summary[idx,3]
idx <- min(intersect(which(spca_summary[,2] >= 500), which(spca_summary[,3] >= 0.9*target_var)))
spca_idx <- sort(unlist(apply(res_list[[idx]]$v, 2, function(x){which(x != 0)})))

descend_idx <- which(colnames(dat) %in% res_hvg$HVG.genes)
gene_idx <- sort(unique(c(spca_idx, descend_idx)))
dat <- dat[,gene_idx]

dat_impute <- dat

###########


dat_se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = t(dat_impute)))
tmp <- zinbwave::zinbwave(dat_se, K = 5, maxiter.optimize = 100, normalizedValues = T,
                          commondispersion = F)
save("../experiment/Writeup_revision7_zinbwave.RData")
