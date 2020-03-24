rm(list=ls())
library(SingleCellExperiment)
source("../simulation/factorization_generator.R")

paramMat <- cbind(50, 120, 10,
                  2, 50, 1/250, 1000,
                  80, 120, 600,
                  1/4, 1/4, 1/2)
colnames(paramMat) <- c("n_each", "d_each", "sigma",
                        "k", "max_iter", "modifier", "max_val",
                        "size_1", "size_2", "size_3",
                        "prop_1", "prop_2", "prop_3")

cell_pop <- matrix(c(4,10, 25,100, 60,80, 25,100,
                     40,10, 60,80, 60,80, 100, 25),
                   nrow = 4, ncol = 4, byrow = T)
gene_pop <- matrix(c(20,90, 25,100,
                     90,20, 100,25)/20, nrow = 2, ncol = 4, byrow = T)

vec <- paramMat[1,]
n_each <- vec["n_each"]
d_each <- vec["d_each"]
sigma <- vec["sigma"]
modifier <- vec["modifier"]

res <- generate_natural_mat(cell_pop, gene_pop, n_each, d_each, sigma, modifier)
nat_mat <- res$nat_mat

r_vec <- sample(c(paramMat[1,"size_1"], paramMat[1,"size_2"], paramMat[1,"size_3"]),
                size = ncol(nat_mat),
                prob = c(paramMat[1,"prop_1"], paramMat[1,"prop_2"], paramMat[1,"prop_3"]),
                replace = T)
# r_vec <- rep(600, ncol(nat_mat))

dat <- generator_zinb_nb(nat_mat, r_vec)
obs_mat <- round(dat$dat * 1000/max(dat$dat))

cluster_labels <- rep(1:4, each = vec["n_each"])

# set.seed(10)
# dat_se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = t(obs_mat)))
# tmp <- zinbwave::zinbwave(dat_se, K = 2, maxiter.optimize = 100, normalizedValues = F,
#                           commondispersion = F)
# res_zinb <- SingleCellExperiment::reducedDims(tmp)$zinbwave
# plot(res_zinb[,1], res_zinb[,2], asp = T, pch = 16, col = cluster_labels)
#
#
# tmp <- zinbwave::zinbwave(dat_se, K = 2, maxiter.optimize = 100, normalizedValues = F,
#                           commondispersion = T)
# res_zinb <- SingleCellExperiment::reducedDims(tmp)$zinbwave
# plot(res_zinb[,1], res_zinb[,2], asp = T, pch = 16, col = cluster_labels)
#
# dist_mat_truth <- as.matrix(stats::dist(res$cell_mat))
# dist_mat_est <- as.matrix(stats::dist(res_zinb))
# mean(sapply(1:nrow(dist_mat_est), function(i){
#   cor(dist_mat_truth[i,], dist_mat_est[i,], method = "kendall")
# }))

###############################

custom.settings = umap.defaults
custom.settings$n_neighbors = 75
custom.settings$min_dist = 0.1
custom.settings$init = "random"

paramMat <- as.matrix(expand.grid(c(2,3,4,5,10,15,20,30,50,75),
                                  c(1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.2, 0.3, 0.5, 0.75, 0.9)))
colnames(paramMat) <- c("n_neighbors", "min_dist")

for(i in 1:nrow(paramMat)){
  custom.settings <- umap.defaults
  custom.settings$n_neighbors <- paramMat[i, "n_neighbors"]
  custom.settings$min_dist <- paramMat[i, "min_dist"]
  custom.settings$init = "random"

  set.seed(10)
  zz <- umap::umap(obs_mat, config = custom.settings)

  png(filename = paste0("../../esvd_results/figure/experiment/Revision_writeup5_umap_", i, ".png"),
      height = 1000, width = 1000, res = 300,
      units = "px")
  plot(zz$layout[,1], zz$layout[,2], asp = T, pch = 16, col = cluster_labels)
  graphics.off()
}

perplexity_vec <- 2:50

for(perplexity in perplexity_vec){
  set.seed(10)
  tmp <- Rtsne::Rtsne(obs_mat, perplexity = perplexity)
  res_tsne <- tmp$Y

  png(filename = paste0("../../esvd_results/figure/experiment/Revision_writeup5_tsne_", perplexity, ".png"),
      height = 1000, width = 1000, res = 300,
      units = "px")
  plot(res_tsne[,1], res_tsne[,2], asp = T, pch = 16, col = cluster_labels)
  graphics.off()
}

##################################

scalar <- 50
init <- eSVD::initialization(obs_mat, family = "neg_binom", k = 3, max_val = 2000,
                             scalar = scalar)
fit <- eSVD::fit_factorization(obs_mat, u_mat = init$u_mat, v_mat = init$v_mat,
                               family = "neg_binom", scalar = scalar,
                               max_iter = 50, max_val = 2000,
                               return_path = F, cores = NA,
                               verbose = T)
plot(fit$u_mat[,1], fit$u_mat[,2], asp = T, pch = 16, col = cluster_labels)

dist_mat_truth <- as.matrix(stats::dist(res$cell_mat))
dist_mat_est <- as.matrix(stats::dist(fit$u_mat[,1:2]))
mean(sapply(1:nrow(dist_mat_est), function(i){
  cor(dist_mat_truth[i,], dist_mat_est[i,], method = "kendall")
}))

nat_mat <- fit$u_mat %*% t(fit$v_mat)
mean_mat <- compute_mean(nat_mat, family = "neg_binom", scalar = scalar)
plot(mean_mat, obs_mat, asp = T, pch = 16, col = rgb(0.5, 0.5, 0.5, 0.2))
lines(c(-1e5,1e5), c(-1e5,1e5), col = "red", lwd = 2, lty = 2)

##################################

set.seed(10)
tmp <- pCMF::pCMF(obs_mat, K = 2, sparsity = F, verbose = F)
res_pcmf <- tmp$factor$U

plot(res_pcmf[,1], res_pcmf[,2], asp = T,  pch = 16, col = cluster_labels)

###########

plot(zz$fit$u_mat[,1], zz$fit$u_mat[,2], pch = 16, col = cluster_labels, asp = T)
