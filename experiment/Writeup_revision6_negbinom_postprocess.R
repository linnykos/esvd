rm(list=ls())
load("../results/factorization_results_negbinom_esvd.RData")
for(i in 1:length(res)){
  len_vec <- sapply(res[[i]], length)
  res[[i]] <- res[[i]][which(len_vec > 1)]
}
for(i in 1:length(res[[1]])){
  res[[1]][[i]]$fit$fit <- res[[1]][[i]]$fit$fit$u_mat
}
res_tmp <- res

load("../results/factorization_results_negbinom_rest.RData")
for(i in 1:length(res)){
  len_vec <- sapply(res[[i]], length)
  res[[i]] <- res[[i]][which(len_vec > 1)]
}
res[[6]] <- res[[5]]
res[3:5] <- res[2:4]
res[2] <- res_tmp[1]

res[c(1,2)] <- res[c(2,1)]

################

# cluster_labels <- rep(1:4, each = 50)
# j <- 3
# i <- 1
# par(mfrow = c(1,2))
# plot(res[[j]][[i]]$fit$fit[,1], res[[j]][[i]]$fit$fit[,2], asp = T,
#      col = cluster_labels, pch = 16, main = "Estimated")
# plot(res[[j]][[i]]$truth[,1], res[[j]][[i]]$truth[,2], asp = T,
#      col = cluster_labels, pch = 16, main = "Truth")
#
# j <- 2
# quality_vec <- sapply(1:length(res[[j]]), function(i){
#   dist_mat_truth <- as.matrix(stats::dist(res[[j]][[i]]$truth))
#   dist_mat_est <- as.matrix(stats::dist(res[[j]][[i]]$fit$fit[,1:2]))
#
#   mean(sapply(1:nrow(dist_mat_est), function(i){
#     cor(dist_mat_truth[i,], dist_mat_est[i,], method = "kendall")
#   }))
# })
# quantile(quality_vec)

####################

trials <- min(sapply(res, length))
res_mat <- matrix(NA, length(res), trials)
for(i in 1:trials){
  if(i %% floor(trials/10) == 0) cat('*')

  for(j in 1:6){
    dist_mat_truth <- as.matrix(stats::dist(res[[j]][[i]]$truth))
    dist_mat_est <- as.matrix(stats::dist(res[[j]][[i]]$fit$fit[,1:2]))

    res_mat[j,i] <- mean(sapply(1:nrow(dist_mat_est), function(x){
      cor(dist_mat_truth[x,], dist_mat_est[x,], method = "kendall")
    }))
  }
}

color_func <- function(alpha = 0.2){
  c(rgb(240/255, 228/255, 66/255, alpha), #yellow
    rgb(86/255, 180/255, 233/255, alpha), #skyblue
    rgb(0/255, 158/255, 115/255, alpha), #bluish green
    rgb(0/255, 114/255, 178/255, alpha), #blue
    rgb(230/255, 159/255, 0/255, alpha), #orange
    rgb(150/255, 150/255, 150/255, alpha))
}

# start of intensive plotting function
den_list <- lapply(1:nrow(res_mat), function(i){
  density(res_mat[i,])
})

#max_val <- max(sapply(den_list, function(x){max(x$y)}))
scaling_factor <- quantile(sapply(den_list, function(x){max(x$y)}), probs = 0.3)

col_vec <- color_func(1)[c(5,2,3,1,4,6)]
text_vec <- c("eSVD", "SVD", "ZINB-WaVE", "pCMF", "(Oracle)\nUMAP", "(Oracle)\nt-SNE")
max_height <- 3

png(paste0("../../esvd_results/figure/experiment/factorization_negbinom_density.png"),
    height = 1800, width = 1000, res = 300, units = "px")
par(mar = c(4,0.5,4,0.5))
plot(NA, xlim = c(-0.3, 1), ylim = c(0, 6.2), ylab = "",
     yaxt = "n", bty = "n", xaxt = "n", xlab = "Kendall's tau",
     main = paste0("Relative embedding correlation"))
axis(side = 1, at = seq(0,1,length.out = 6))
for(i in 1:nrow(res_mat)){
  lines(c(0,1), rep(nrow(res_mat) - i, 2))

  y_vec <- (c(0, den_list[[i]]$y, 0 , 0))/scaling_factor
  if(max(y_vec) > max_height) y_vec <- y_vec*max_height/max(y_vec)
  polygon(x = c(den_list[[i]]$x[1], den_list[[i]]$x, den_list[[i]]$x[length(den_list[[i]]$x)], den_list[[i]]$x[1]),
          y = y_vec + nrow(res_mat) - i,
          col = col_vec[i])

  med <- median(res_mat[i,])
  lines(rep(med, 2), y = c(nrow(res_mat) - i, 0), lwd = 1, lty = 2)
  points(med, y = nrow(res_mat) - i, col = "black", pch = 16, cex = 2)
  points(med, y = nrow(res_mat) - i, col = col_vec[i], pch = 16, cex = 1.5)
}
text(x = rep(-0.15,6), y = seq(5.35, 0.35, by=-1), labels = text_vec)
graphics.off()

##########################


col_func2 <- function(alpha){
  c( rgb(86/255, 180/255, 233/255, alpha), #skyblue
     rgb(240/255, 228/255, 66/255, alpha), #yellow
     rgb(0/255, 158/255, 115/255, alpha), #bluish green
     rgb(230/255, 159/255, 0/255,alpha)) #orange
}
col_vec <- col_func2(1)

png(paste0("../../esvd_results/figure/experiment/factorization_negbinom_embedding.png"),
    height = 1500, width = 1500, res = 300, units = "px")
text_vec <- c("eSVD", "SVD", "ZINB-WaVE", "pCMF", "(Oracle) UMAP", "(Oracle) t-SNE")
par(mfrow = c(2,3), mar = c(1, 1, 1.5, 1))
for(i in 1:6){
  idx <- which.min(abs(res_mat[i,] - median(res_mat[i,])))

  plot(res[[i]][[idx]]$fit$fit[,1], res[[i]][[idx]]$fit$fit[,2],
       asp = T, pch = 16, col = col_vec[rep(1:4, each = paramMat[1,"n_each"])],
       xlab = "Latent dim. 1", ylab = "Latent dim. 2",
       main = paste0(text_vec[i], ": (", round( median(res_mat[i,]), 2), ")"),
       xaxt = "n", yaxt = "n")
}
graphics.off()
