rm(list=ls())
load("../results/factorization_results_curved_gaussian_esvd.RData")
for(i in 1:length(res)){
  len_vec <- sapply(res[[i]], length)
  res[[i]] <- res[[i]][which(len_vec > 1)]
}
for(i in 1:length(res[[1]])){
  res[[1]][[i]]$fit$fit <- res[[1]][[i]]$fit$fit$u_mat
}
res_tmp <- res

load("../results/factorization_results_curved_gaussian_rest.RData")
for(i in 1:length(res)){
  len_vec <- sapply(res[[i]], length)
  res[[i]] <- res[[i]][which(len_vec > 1)]
}
res[[length(res)+1]] <- res_tmp[[1]]

current_ord <- c("SVD", "ZINB-WaVE", "pCMF", "(Oracle)\nUMAP", "(Oracle)\nt-SNE",
                 "Isomap", "ICA", "NMF", "Diffusion\nmap", "eSVD (Curved\nGaussian)")
desired_ord <- c("eSVD (Curved\nGaussian)",
                 "ZINB-WaVE", "pCMF",
                 "SVD", "NMF", "ICA",
                 "(Oracle)\nUMAP", "(Oracle)\nt-SNE", "Isomap", "Diffusion\nmap")

res_tmp <- res
for(i in 1:length(desired_ord)){
  idx <- which(current_ord == desired_ord[i])
  res[[i]] <- res_tmp[[idx]]
}

trials <- min(sapply(res, length))
res_mat <- matrix(NA, length(res), trials)
for(i in 1:trials){
  if(i %% floor(trials/10) == 0) cat('*')

  for(j in 1:length(res)){
    dist_mat_truth <- as.matrix(stats::dist(res[[j]][[i]]$truth))
    dist_mat_est <- as.matrix(stats::dist(res[[j]][[i]]$fit$fit[,1:2]))

    res_mat[j,i] <- mean(sapply(1:nrow(dist_mat_est), function(x){
      cor(dist_mat_truth[x,], dist_mat_est[x,], method = "kendall")
    }))
  }
}

####################

color_func <- function(alpha = 0.2){
  c(rgb(240/255, 228/255, 66/255, alpha), #yellow 1
    rgb(86/255, 180/255, 233/255, alpha), #skyblue 2
    rgb(0/255, 158/255, 115/255, alpha), #bluish green 3
    rgb(0/255, 114/255, 178/255, alpha), #blue 4
    rgb(230/255, 159/255, 0/255, alpha), #orange 5
    rgb(150/255, 150/255, 150/255, alpha), # gray 6
    rgb(204/255, 121/255, 167/255, alpha), # pink 7
    rgb(213/255, 94/255, 0/255, alpha), #dark orange 8
    rgb(67/255, 24/255, 97/255, alpha), #dark purple 9
    rgb(25/255, 67/255, 32/255, alpha) #dark green 10
  )
}

# start of intensive plotting function
den_list <- lapply(1:nrow(res_mat), function(i){
  density(res_mat[i,])
})

#max_val <- max(sapply(den_list, function(x){max(x$y)}))
scaling_factor <- quantile(sapply(den_list, function(x){max(x$y)}), probs = 0.35)

col_vec <- color_func(1)[c(5,2,3,1,4,8,6,10,7,9)]
text_vec <- desired_ord
max_height <- 1
y_spacing <- .5

png(paste0("../../esvd_results/figure/simulation/factorization_curved_gaussian_density.png"),
    height = 2500, width = 1000, res = 300, units = "px")
par(mar = c(4,0.5,4,0.5))
plot(NA, xlim = c(-0.3, 1), ylim = c(0, y_spacing*nrow(res_mat)+0.2), ylab = "",
     yaxt = "n", bty = "n", xaxt = "n", xlab = "Kendall's tau",
     main = paste0("Relative embedding correlation:\nCurved Gaussian generative model"),
     cex.main = 1.1)
axis(side = 1, at = seq(0,1,length.out = 6))
for(i in 1:nrow(res_mat)){
  lines(c(0,1), rep((nrow(res_mat) - i)*y_spacing, 2))

  y_vec <- (c(0, den_list[[i]]$y, 0 , 0))/scaling_factor
  if(max(y_vec) > max_height) y_vec <- y_vec*max_height/max(y_vec)
  polygon(x = c(den_list[[i]]$x[1], den_list[[i]]$x, den_list[[i]]$x[length(den_list[[i]]$x)], den_list[[i]]$x[1]),
          y = y_vec + (nrow(res_mat) - i)*y_spacing,
          col = col_vec[i])

  med <- median(res_mat[i,])
  lines(rep(med, 2), y = c((nrow(res_mat) - i)*y_spacing, 0), lwd = 1, lty = 2)
  points(med, y = (nrow(res_mat) - i)*y_spacing, col = "black", pch = 16, cex = 2)
  points(med, y = (nrow(res_mat) - i)*y_spacing, col = col_vec[i], pch = 16, cex = 1.5)
}

text(x = rep(-0.1, nrow(res_mat)), y = seq((nrow(res_mat)-1+.35)*y_spacing, 0.35*y_spacing, by = -y_spacing), labels = text_vec,
     cex = 0.9)
graphics.off()

##########################


col_func2 <- function(alpha){
  c( rgb(86/255, 180/255, 233/255, alpha), #skyblue
     rgb(240/255, 228/255, 66/255, alpha), #yellow
     rgb(0/255, 158/255, 115/255, alpha), #bluish green
     rgb(230/255, 159/255, 0/255,alpha)) #orange
}
col_vec <- col_func2(1)

set.seed(10)
shuff_idx <- sample(1:nrow(res[[1]][[1]]$fit$fit))

png(paste0("../../esvd_results/figure/simulation/factorization_curved_gaussian_embedding.png"),
    height = 1500, width = 1800, res = 300, units = "px")
text_vec <- c("eSVD", "ZINB-WaVE", "pCMF", "SVD", "NMF", "ICA",
              "UMAP", "t-SNE", "Isomap", "Diff. Map")
par(mfrow = c(2,5), mar = c(0.5, 0.5, 1.5, 0.25))
for(i in 1:nrow(res_mat)){
  idx <- which.min(abs(res_mat[i,] - median(res_mat[i,])))

  tmp_mat <- res[[i]][[idx]]$fit$fit[,1:2]

  # rotate the non-eSVD factorizations to fit with the plot better
  if(i != 1){
    tmp_mat <- scale(tmp_mat, center = T, scale = F)

    idx_esvd <-  which.min(abs(res_mat[1,] - median(res_mat[1,])))
    tmp_esvd <- res[[1]][[idx_esvd]]$fit$fit[,1:2]
    tmp_esvd <- scale(tmp_esvd, center = T, scale = T)

    tmp_svd <- svd(t(tmp_mat) %*% tmp_esvd)
    rot_mat <- tmp_svd$u[,1:2] %*% t(tmp_svd$v[,1:2])
    tmp_mat <- tmp_mat %*% rot_mat
  }

  xlim <- range(tmp_mat[,1]); ylim <- range(tmp_mat[,2])

  # construct the 5-column grid points
  x_seq <- seq(xlim[1], xlim[2], length.out = 5)
  diff_val <- abs(x_seq[2] - x_seq[1])
  y_seq <- c(-10:10)*diff_val + mean(ylim)

  plot(NA, xlim = xlim, ylim = ylim,
       asp = T,
       xlab = "Latent dim. 1", ylab = "Latent dim. 2",
       main = paste0(text_vec[i], ": (", round( median(res_mat[i,]), 2), ")"),
       xaxt = "n", yaxt = "n")

  # plot grid
  for(i in 1:length(x_seq)){
    lines(rep(x_seq[i], 2), c(-1e5,1e5), lwd = 0.5, lty = 3, col = "gray")
  }
  for(i in 1:length(y_seq)){
    lines(c(-1e5,1e5), rep(y_seq[i], 2), lwd = 0.5, lty = 3, col = "gray")
  }

  for(i in 1:nrow(tmp_mat)){
    points(tmp_mat[shuff_idx[i],1], tmp_mat[shuff_idx[i],2],
           pch = 16, cex = 1.2)
    points(tmp_mat[shuff_idx[i],1], tmp_mat[shuff_idx[i],2],
           pch = 16, col = col_vec[rep(1:4, each = paramMat[1,"n_each"])][shuff_idx[i]])
  }

}
graphics.off()
