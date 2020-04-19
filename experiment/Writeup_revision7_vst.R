rm(list=ls())
load("../results/step5_trajectory_log2.RData")
gene_mean <- apply(dat_impute, 2, mean)
gene_var <- apply(dat_impute, 2, sd)

par(mfrow = c(1,2))
plot(gene_mean, gene_var, asp = T,
     xlab = "Mean (per gene)", ylab = "Sd (per gene)", pch = 16)
lines(c(-1e4,1e4), c(-1e4,1e4), col = "red", lty = 2)
plot(log(gene_mean), log(gene_var), asp = T,
     xlab = "Log of mean (per gene)", ylab = "Log of sd (per gene)", pch = 16)
lines(c(-1e4,1e4), c(-1e4,1e4), col = "red", lty = 2)
