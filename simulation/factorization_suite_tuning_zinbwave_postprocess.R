rm(list=ls())
load("../results/factorization_results_tuning_zinbwave.RData")

# # see if there's an empirical difference
# x <- 1; i <- 3
# plot(res[[i]][[x]]$fit[[1]]$u_mat[,1], res[[i]][[x]]$fit[[1]]$u_mat[,2], asp = T, pch = 16, col = rep(1:4, each = 50))

loglik_list <- lapply(1:trials, function(x){
  if(x %% floor(trials/10) == 0) cat('*')

  nat_mat_list_list <- lapply(1:nrow(paramMat), function(i){
    lapply(1:3, function(j){
      u_mat <- res[[i]][[x]]$fit[[j]]$u_mat
      v_mat <- res[[i]][[x]]$fit[[j]]$v_mat
      u_mat %*% t(v_mat)
    })
  })

  loglik_vec <- sapply(1:nrow(paramMat), function(i){
    tmp_vec <- sapply(1:3, function(j){
      idx <- res[[1]][[x]]$missing_idx[[j]]
      dat_vec <- res[[1]][[x]]$dat[idx]
      nat_vec <- nat_mat_list_list[[i]][[j]][idx]
      r_val <- paramMat[i,"r_val"]

      sum(-r_val*log(1-exp(nat_vec)) - dat_vec*nat_vec)/length(idx)
    })

    mean(tmp_vec)
  })
})

table(sapply(loglik_list, which.min))
loglik_mat <- do.call(cbind, loglik_list)

# evaluate the accuracy of the fits
accuracy_mat <- matrix(NA, nrow(paramMat), trials)

for(i in 1:trials){
  if(i %% floor(trials/10) == 0) cat('*')

  for(j in 1:nrow(paramMat)){
    accuracy_mat[j,i] <- mean(sapply(1:3, function(k){
      dist_mat_truth <- as.matrix(stats::dist(res[[j]][[i]]$true_u_mat))
      dist_mat_est <- as.matrix(stats::dist(res[[j]][[i]]$fit[[k]]$u_mat[,1:2]))

      mean(sapply(1:nrow(dist_mat_est), function(x){
        cor(dist_mat_truth[x,], dist_mat_est[x,], method = "kendall")
      }))
    }))
  }
}

table(apply(accuracy_mat, 2, which.max))


png(paste0("../../esvd_results/figure/simulation/factorization_tuning.png"),
    height = 1000, width = 2500, res = 300, units = "px")

par(mfrow = c(1,4), mar = c(4,4,4,0.2))
plot(as.numeric(loglik_mat), as.numeric(accuracy_mat), pch = 16, xlab = "Negative log-likelihood:\nTesting set",
     ylab = "Relative embedding correlation", main = "Performance across all\ntrials & tuning parameters",
     cex.main = 1.1)
idx <- apply(loglik_mat, 2, which.min)
for(i in 1:length(idx)){
  col <- ifelse(idx[i] == 2, rgb(86/255, 180/255, 233/255), rgb(230/255, 159/255, 0/255))
  points(loglik_mat[idx[i],i], accuracy_mat[idx[i],i], pch = 16, cex = 1.4)
  points(loglik_mat[idx[i],i], accuracy_mat[idx[i],i], pch = 16, col = col, cex = 1.2)
}

legend("bottomright", c("Selected, k=3", "Selected, k=2", "Others"),
       fill = c(rgb(86/255, 180/255, 233/255), rgb(230/255, 159/255, 0/255), "black"),
       cex = 0.8)

col_func2 <- function(alpha){
  c( rgb(86/255, 180/255, 233/255, alpha), #skyblue
     rgb(240/255, 228/255, 66/255, alpha), #yellow
     rgb(0/255, 158/255, 115/255, alpha), #bluish green
     rgb(230/255, 159/255, 0/255,alpha)) #orange
}
col_vec <- col_func2(1)

par(mar = c(4,1.5,4,0.2))
for(i in 1:3){
  idx <- which.min(abs(accuracy_mat[i,] - median(accuracy_mat[i,])))

  tmp_mat <- res[[i]][[idx]]$fit[[1]]$u_mat[,1:2]

  # rotate the non-eSVD factorizations to fit with the plot better
  if(i != 2){
    tmp_mat <- scale(tmp_mat, center = T, scale = F)

    idx_k2 <-  which.min(abs(accuracy_mat[2,] - median(accuracy_mat[2,])))
    tmp_esvd <- res[[2]][[idx_k2]]$fit[[1]]$u_mat[,1:2]
    tmp_esvd <- scale(tmp_esvd, center = T, scale = T)

    tmp_svd <- svd(t(tmp_mat) %*% tmp_esvd)
    rot_mat <- tmp_svd$u[,1:2] %*% t(tmp_svd$v[,1:2])
    tmp_mat <- tmp_mat %*% rot_mat
  }

  xlim <- range(tmp_mat[,1]); ylim <- range(tmp_mat[,2])

  # construct the 5-column grid points
  if(diff(xlim) > diff(ylim)){
    x_seq <- seq(xlim[1], xlim[2], length.out = 9)
    diff_val <- abs(x_seq[2] - x_seq[1])
    y_seq <- c(-15:15)*diff_val + mean(ylim)
  } else {
    y_seq <- seq(ylim[1], ylim[2], length.out = 9)
    diff_val <- abs(y_seq[2] - y_seq[1])
    x_seq <- c(-15:15)*diff_val + mean(xlim)
  }

  plot(NA, xlim = xlim, ylim = ylim,
       asp = T,
       xlab = "Latent dim. 1", ylab = "",
       main = paste0("k = ", paramMat[i,"k"], ": Relative embedding\ncorrelation = ", round( median(accuracy_mat[i,]), 2)),
       xaxt = "n", yaxt = "n", cex.main = 1.1)
  title(ylab="Latent dim. 2", line = .25)

  # plot grid
  for(i in 1:length(x_seq)){
    lines(rep(x_seq[i], 2), c(-1e5,1e5), lwd = 0.5, lty = 3, col = "gray")
  }
  for(i in 1:length(y_seq)){
    lines(c(-1e5,1e5), rep(y_seq[i], 2), lwd = 0.5, lty = 3, col = "gray")
  }

  for(i in 1:nrow(tmp_mat)){
    points(tmp_mat[i,1], tmp_mat[i,2], pch = 16, cex = 1.2)
    points(tmp_mat[i,1], tmp_mat[i,2],  pch = 16, col = col_vec[rep(1:4, each = paramMat[1,"n_each"])][i])
  }

}
graphics.off()

