rm(list=ls())
load("../results/factorization_results.RData")

cluster_labels <- c(1:4)[rep(1:4, each = paramMat[1,"n"])]
.b_estimate <- function(mat, cluster_labels){
  uniq_val <- sort(unique(cluster_labels))
  mean(sapply(uniq_val, function(x){
    idx <- which(cluster_labels == x)
    stats::princomp(mat[idx,])$sdev[1]
  }))/4
}

.compare_two_lineages <- function(target_lineage, source_lineage, cluster_labels){
  len1 <- length(target_lineage)
  len2 <- length(source_lineage)

  sapply(1:len1, function(x){
    target_lambda <- target_lineage[[x]]$lambda_long
    idx <- which(target_lambda != 0)
    target_lambda <- target_lambda[idx]

    max(sapply(1:len2, function(y){
      source_lambda <- source_lineage[[y]]$lambda_long
      source_lambda <- source_lambda[idx]
      abs(stats::cor(target_lambda, source_lambda))
    }))
  })
}


# temporarily, use the ZINF-WAVE model
res_zinb_list <- vector("list", length = trials)
for(i in 1:trials){
  print(i)
  dat_se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = t(res[[1]][[1]]$dat)))
  tmp <- zinbwave::zinbwave(dat_se, K = 4, maxiter.optimize = 100)
  res_zinb_list[[i]] <- tmp@reducedDims$zinbwave
}
# save.image("../results/factorization_results2.RData")


b_est <- .b_estimate(res[[1]][[1]]$cell_mat[,1:2], cluster_labels)
fixed_lineage <- singlecell::slingshot(res[[1]][[1]]$cell_mat[,1:2], cluster_labels, 1, knn = NA,
                                       b = b_est, remove_outlier = F)$lineages

kendall_list <- vector("list", length(res[[1]]))
for(x in 1:length(kendall_list)){
  print(x)
  set.seed(x)
  b_est <- .b_estimate(res[[1]][[x]]$cell_mat[,1:2], cluster_labels)
  true_lineage <- .get_curves(res[[1]][[x]]$cell_mat[,1:2], cluster_labels, fixed_lineage, shrink = 1,
                              thresh =  0.001, max_iter = 15, b = b_est)

  # plot(res[[1]][[x]]$cell_mat[,1], res[[1]][[x]]$cell_mat[,2], asp = T, pch = 16, col = rgb(0,0,0,0.1))
  # for(i in 1:length(true_lineage)){
  #   ord <- true_lineage[[i]]$ord
  #   lines(true_lineage[[i]]$s[ord, 1], true_lineage[[i]]$s[ord, 2], lwd = 3,
  #         col = "black")
  # }

  b_est <- .b_estimate(res[[1]][[x]]$res_our$u_mat[,1:2], cluster_labels)
  our_lineage <- .get_curves(res[[1]][[x]]$res_our$u_mat[,1:2], cluster_labels, fixed_lineage, shrink = 1,
                             thresh =  0.001, max_iter = 15, b = b_est)

  # plot(res[[1]][[x]]$res_our$u_mat[,1], res[[1]][[x]]$res_our$u_mat[,2], asp = T, pch = 16, col = cluster_labels)
  # for(i in 1:length(our_lineage)){
  #   ord <- our_lineage[[i]]$ord
  #   lines(our_lineage[[i]]$s[ord, 1], our_lineage[[i]]$s[ord, 2], lwd = 3,
  #         col = "black")
  # }

  b_est <- .b_estimate(res[[1]][[x]]$res_svd[,1:2], cluster_labels)
  svd_lineage <- .get_curves(res[[1]][[x]]$res_svd[,1:2], cluster_labels, fixed_lineage, shrink = 1,
                             thresh =  0.001, max_iter = 15, b = b_est)

  # plot(res[[1]][[x]]$res_svd[,1], res[[1]][[x]]$res_svd[,2], asp = T, pch = 16, col = cluster_labels)
  # for(i in 1:length(svd_lineage)){
  #   ord <- svd_lineage[[i]]$ord
  #   lines(svd_lineage[[i]]$s[ord, 1], svd_lineage[[i]]$s[ord, 2], lwd = 3,
  #         col = "black")
  # }

  b_est <- .b_estimate(res[[1]][[x]]$res_ica[,1:2], cluster_labels)
  ica_lineage <- singlecell::slingshot(res[[1]][[x]]$res_ica[,1:2],
                                       cluster_labels, 1, knn = NA,
                                       b = b_est, remove_outlier = F)$curves

  # estimate the TSNE embedding too, since I forgot to do this
  res_tsne <-  Rtsne::Rtsne(res[[1]][[x]]$dat_impute, dims = 2, perplexity=30, verbose=F, max_iter = 500)
  b_est <- .b_estimate(res_tsne$Y, cluster_labels)
  tsne_lineage <- singlecell::slingshot(res_tsne$Y,
                                        cluster_labels, 1, knn = NA,
                                        b = b_est, remove_outlier = F)$curves

  b_est <- .b_estimate(res_zinb_list[[x]], cluster_labels)
  zinf_lineage <- singlecell::slingshot(res_zinb_list[[x]][,1:2], cluster_labels, 1, knn = NA,
                                        b = b_est, remove_outlier = F)$curves

  res_mat <- matrix(c(.compare_two_lineages(true_lineage, our_lineage),
                      .compare_two_lineages(true_lineage, svd_lineage),
                      .compare_two_lineages(true_lineage, ica_lineage),
                      .compare_two_lineages(true_lineage, tsne_lineage),
                      .compare_two_lineages(true_lineage, zinf_lineage)),
                    ncol = 2, byrow = T)
  rownames(res_mat) <- c("Our", "SVD", "ICA", "tSNE", "ZINF")
  kendall_list[[x]] <- res_mat
}

our_vec <- unlist(lapply(kendall_list, function(x){x["Our",]}))
svd_vec <- unlist(lapply(kendall_list, function(x){x["SVD",]}))
ica_vec <- unlist(lapply(kendall_list, function(x){x["ICA",]}))
tsne_vec <- unlist(lapply(kendall_list, function(x){x["tSNE",]}))
zinf_vec <- unlist(lapply(kendall_list, function(x){x["ZINF",]}))

# hist(our_vec, breaks = 25, col = "gray", xlim = c(0,1))
# hist(svd_vec, breaks = 25, col = "gray", xlim = c(0,1))
# hist(ica_vec, breaks = 25, col = "gray", xlim = c(0,1))
# hist(tsne_vec, breaks = 25, col = "gray", xlim = c(0,1))

d_list <- list(density(our_vec), density(svd_vec), density(ica_vec), density(tsne_vec),
               density(zinf_vec))
d_list <- lapply(d_list, function(x){
  idx <- which(x$x<=1)
  x$x <- c(x$x[idx],1)
  x$y <- c(x$y[idx],0)
  x
})
d_list[[1]]$x <- c(d_list[[1]]$x, 1); d_list[[1]]$y <- c(d_list[[1]]$y, 0)
ylim <- range(unlist(lapply(d_list, function(x){x$y})))
red <- rgb(205,40,54,maxColorValue=255)

png("../figure/simulation/factorization_density_2.png", height = 800, width = 2500, res = 300, units = "px")
par(mfrow = c(1,5), mar = c(4,4,4,0.5))
plot(d_list[[1]], main=paste0("Our method (", round(median(our_vec), 2), ")"), ylab = "Density", xlab = "Pearson correlation",
     xlim = c(0,1), cex.lab = 1.25)
polygon(d_list[[1]], col="gray", border="black")
lines(rep(median(our_vec), 2), c(0, 1e5), lwd = 4, col = red, lty = 2)

plot(d_list[[2]],main=paste0("SVD (", round(median(svd_vec), 2), ")"), ylab = "Density", xlab = "Pearson correlation",
     xlim = c(0,1), cex.lab = 1.25)
polygon(d_list[[2]], col="gray", border="black")
lines(rep(median(svd_vec), 2), c(0, 1e5), lwd = 4, col = red, lty = 2)

plot(d_list[[3]], main=paste0("ICA (", round(median(ica_vec), 2), ")"), ylab = "Density", xlab = "Pearson correlation",
     xlim = c(0,1), cex.lab = 1.25)
polygon(d_list[[3]], col="gray", border="black")
lines(rep(median(ica_vec), 2), c(0, 1e5), lwd = 4, col = red, lty = 2)

plot(d_list[[4]], main=paste0("t-SNE (", round(median(tsne_vec), 2), ")"), ylab = "Density", xlab = "Pearson correlation",
     xlim = c(0,1), cex.lab = 1.25)
polygon(d_list[[4]], col="gray", border="black")
lines(rep(median(tsne_vec), 2), c(0, 1e5), lwd = 4, col = red, lty = 2)

plot(d_list[[5]], main=paste0("ZINF-WaVE (", round(median(zinf_vec), 2), ")"), ylab = "Density", xlab = "Pearson correlation",
     xlim = c(0,1), cex.lab = 1.25)
polygon(d_list[[5]], col="gray", border="black")
lines(rep(median(zinf_vec), 2), c(0, 1e5), lwd = 4, col = red, lty = 2)
graphics.off()

#####################

# make the example plot first

col_vec <- c(rgb(205,40,54,maxColorValue=255), #red
             rgb(180,200,255,maxColorValue=255), #blue
             rgb(100,100,200,maxColorValue=255), #purple
             rgb(149,219,144,maxColorValue=255)) #green

x <- 4
set.seed(10)
png("../figure/simulation/factorization_illustration.png", height = 600, width = 1900, res = 300, units = "px")
set.seed(10)
par(mfrow = c(1,5), mar = c(0.5, 0.5, 3, 0.5))
# x <- ceiling(which.min(abs(svd_vec - 0.8))/2)
# b_est <- .b_estimate(res[[1]][[x]]$cell_mat[,1:2], cluster_labels)
# true_lineage <- .get_curves(res[[1]][[x]]$cell_mat[,1:2], cluster_labels, fixed_lineage, shrink = 1,
#                             thresh =  0.001, max_iter = 15, b = b_est)
# plot(res[[1]][[x]]$cell_mat[,1], res[[1]][[x]]$cell_mat[,2], asp = T,
#      pch = 16, col = col_vec[rep(1:4, each = paramMat[1,"n"])], xlab = "", ylab = "",
#      xaxt = "n", yaxt = "n", main = "Truth")
# for(i in 1:length(true_lineage)){
#   ord <- true_lineage[[i]]$ord
#   lines(true_lineage[[i]]$s[ord, 1], true_lineage[[i]]$s[ord, 2], lwd = 2,
#         col = "black")
# }

b_est <- .b_estimate(res[[1]][[x]]$res_our$u_mat[,1:2], cluster_labels)
our_lineage <- .get_curves(res[[1]][[x]]$res_our$u_mat[,1:2], cluster_labels, fixed_lineage, shrink = 1,
                           thresh =  0.001, max_iter = 15, b = b_est)
plot(res[[1]][[x]]$res_our$u_mat[,1], res[[1]][[x]]$res_our$u_mat[,2], asp = T,
     pch = 16, col = col_vec[rep(1:4, each = paramMat[1,"n"])], xlab = "", ylab = "",
     xaxt = "n", yaxt = "n", main = paste0("Our method (", round(mean(.compare_two_lineages(true_lineage, our_lineage)), 2),
                                           ")"))
for(i in 1:length(our_lineage)){
  ord <- our_lineage[[i]]$ord
  lines(our_lineage[[i]]$s[ord, 1], our_lineage[[i]]$s[ord, 2], lwd = 2,
        col = "black")
}

b_est <- .b_estimate(res[[1]][[x]]$res_svd[,1:2], cluster_labels)*10
svd_lineage <- .get_curves(res[[1]][[x]]$res_svd[,1:2], cluster_labels, fixed_lineage, shrink = 1,
                           thresh =  0.001, max_iter = 15, b = b_est)
plot(-res[[1]][[x]]$res_svd[,1], -res[[1]][[x]]$res_svd[,2], asp = T,
     pch = 16, col = col_vec[rep(1:4, each = paramMat[1,"n"])], xlab = "", ylab = "",
     xaxt = "n", yaxt = "n", main = paste0("SVD (", round(mean(.compare_two_lineages(true_lineage, svd_lineage)), 2),
                                           ")"))
for(i in 1:length(svd_lineage)){
  ord <- svd_lineage[[i]]$ord
  lines(-svd_lineage[[i]]$s[ord, 1], -svd_lineage[[i]]$s[ord, 2], lwd = 2,
        col = "black")
}

b_est <- .b_estimate(res[[1]][[x]]$res_ica[,1:2], cluster_labels)*10
ica_lineage <- singlecell::slingshot(res[[1]][[x]]$res_ica[,1:2],
                                     cluster_labels, 1, knn = NA,
                                     b = b_est, remove_outlier = F)$curves
plot(res[[1]][[x]]$res_ica[,1], res[[1]][[x]]$res_ica[,2], asp = T,
     pch = 16, col = col_vec[rep(1:4, each = paramMat[1,"n"])], xlab = "", ylab = "",
     xaxt = "n", yaxt = "n", main = paste0("ICA (", round(mean(.compare_two_lineages(true_lineage, ica_lineage)), 2),
                                           ")"))
for(i in 1:length(ica_lineage)){
  ord <- ica_lineage[[i]]$ord
  lines(ica_lineage[[i]]$s[ord, 1], ica_lineage[[i]]$s[ord, 2], lwd = 2,
        col = "black")
}

res_tsne <-  Rtsne::Rtsne(res[[1]][[x]]$dat_impute, dims = 2, perplexity=5, verbose=F, max_iter = 500)
b_est <- .b_estimate(res_tsne$Y, cluster_labels)*10
tsne_lineage <- singlecell::slingshot(res_tsne$Y,
                                      cluster_labels, 1, knn = NA,
                                      b = b_est, remove_outlier = F)$curves
plot(res_tsne$Y[,2], -res_tsne$Y[,1], asp = T, pch = 16,
     col = col_vec[rep(1:4, each = paramMat[1,"n"])], xlab = "", ylab = "", xaxt = "n",
     yaxt = "n", main = paste0("t-SNE (", round(mean(.compare_two_lineages(true_lineage, tsne_lineage)), 2),
                                      ")"))
for(i in 1:length(tsne_lineage)){
  ord <- tsne_lineage[[i]]$ord
  lines(tsne_lineage[[i]]$s[ord, 2], -tsne_lineage[[i]]$s[ord, 1], lwd = 2,
        col = "black")
}

b_est <- .b_estimate(res_zinb_list[[x]][,1:2], cluster_labels)*10
zinf_lineage <- singlecell::slingshot(res_zinb_list[[x]][,1:2], cluster_labels, 1, knn = NA,
                                      b = b_est, remove_outlier = F)$curves
plot(res_zinb_list[[x]][,1], -res_zinb_list[[x]][,2], asp = T,
     pch = 16, col = col_vec[rep(1:4, each = paramMat[1,"n"])], xlab = "", ylab = "",
     xaxt = "n", yaxt = "n", main = paste0("ZINF-WaVE (", round(mean(.compare_two_lineages(true_lineage, zinf_lineage)), 2),
                                           ")"))
for(i in 1:length(zinf_lineage)){
  ord <- zinf_lineage[[i]]$ord
  lines(zinf_lineage[[i]]$s[ord, 1], -zinf_lineage[[i]]$s[ord, 2], lwd = 2,
        col = "black")
}
graphics.off()

###########

# sparsity vec
sparsity_vec <- sapply(1:length(kendall_list), function(x){
  length(which(res[[1]][[x]]$dat == 0))/prod(dim(res[[1]][[x]]$dat))
})

