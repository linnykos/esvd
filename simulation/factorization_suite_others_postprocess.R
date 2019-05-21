rm(list=ls())
load("../results/factorization_results_others.RData")

col_func <- function(alpha){
  c( rgb(86/255, 180/255, 233/255, alpha), #skyblue
     rgb(240/255, 228/255, 66/255, alpha), #yellow
     rgb(0/255, 158/255, 115/255, alpha), #bluish green
     rgb(230/255, 159/255, 0/255,alpha)) #orange
}
col <- col_func(1)

cluster_labels <- rep(1:4, each = paramMat[1,"n_each"])

type_vec <- c("gaussian", "curved_gaussian", "zinb", "pcmf")
title_vec <- c("Gaussian", "Curved Gaussian", "Negative Binomial", "Poisson")
for(k in 1:length(type_vec)){
  png(paste0("../figure/simulation/factorization_", type_vec[k], "_example.png"),
      height = 1500, width = 2000, res = 300, units = "px")
  label_vec <- c("SVD", "ICA", "t-SNE", "Our", "ZINB-WaVE", "pCMF")
  idx <- 1

  curves_truth <- singlecell::slingshot(res[[k]][[idx]]$dat$truth, cluster_labels,
                                        starting_cluster = 1,
                                        verbose = F)

  par(mfrow = c(2,3), mar = c(4, 4, 4, 0.5))
  for(i in 1:6){
    curves <- singlecell::slingshot(res[[k]][[idx]][[i]], cluster_labels,
                                    starting_cluster = 1, verbose = F,
                                    use_initialization = F)

    ###

    val <- mean(sapply(1:2, function(j){
      idx <- which(curves_truth$curves[[j]]$lambda_long != 0)

      # find which curve is most suitable
      tmp <- sapply(1:length(curves$curves), function(l){
        cor(curves_truth$curves[[j]]$lambda_long[idx],
            curves$curves[[l]]$lambda_long[idx], method = "kendall")
      })

      max(tmp)
    }))

    ###

    plot(res[[k]][[idx]][[i]][,1], res[[k]][[idx]][[i]][,2],
         asp = T, pch = 16, col = col[rep(1:4, each = paramMat[1,"n_each"])],
         xlab = "Latent dim. 1", ylab = "Latent dim. 2",
         main = paste0(label_vec[i], " : (", round(val, 2), ")\n(", title_vec[k], ")"))

    for(j in 1:length(curves$curves)){
      ord <- curves$curves[[j]]$ord
      lines(curves$curves[[j]]$s[ord,1], curves$curves[[j]]$s[ord,2], lwd = 3, col = "white")
      lines(curves$curves[[j]]$s[ord,1], curves$curves[[j]]$s[ord,2], lwd = 2)}


  }

  graphics.off()
}
