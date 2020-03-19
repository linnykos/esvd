rm(list=ls())

load("../results/factorization_esvd.RData")

################################

correct_idx <- c(5, 18, 31, 48)
scalar_vec <- sapply(1:4, function(k){
  i <- correct_idx[k]
  median(sapply(res[[i]], function(y){
    y$fitting_param
  }))
})
scalar_vec

j <- 1
main_vec <- c("Gaussian\n(constant variance)", "Poisson",
              paste0("Negative binomial\n(true size = 50, est. size = ", round(res[[correct_idx[3]]][[j]]$fitting_param, 2), ")"),
              paste0("Curved Gaussian\n(true alpha = 2, est. alpha = ", round(res[[correct_idx[4]]][[j]]$fitting_param, 2), ")"))
family_vec <- c("gaussian", "poisson", "neg_binom", "curved_gaussian")

png(filename = "../../esvd_results/figure/experiment/Revision_writeup4_simulation_testing.png",
    height = 2250, width = 2250, res = 300,
    units = "px")
par(mfrow = c(2,2))
for(k in 1:4){
  i <- correct_idx[k]

  plot_prediction_against_observed(dat = res[[i]][[j]]$dat,
                                   nat_mat = res[[i]][[j]]$fit_u_mat %*% t(res[[i]][[j]]$fit_v_mat),
                                   family = family_vec[k], scalar = res[[i]][[j]]$fitting_param,
                                   main = main_vec[k])
}

graphics.off()

####################################

png(filename = "../../esvd_results/figure/experiment/Revision_writeup4_simulation_fit.png",
    height = 2250, width = 2250, res = 300,
    units = "px")
pos_vec <- c("topleft", "topleft", "topright", "topright")
par(mfrow = c(2,2))
for(k in 1:4){
  i <- correct_idx[k]
  plot(res[[i]][[1]]$fit_u_mat[,1], res[[i]][[1]]$fit_u_mat[,2], col = rep(1:4, each = paramMat[1,"n_each"]), pch = 16,
       xlab = "Latent dimension 1", ylab = "Latent dimension 2", main =
         paste0("Estimated cell latent positions\n(", main_vec2[k], ")"), asp = T)
  legend(pos_vec[k], c("Cell type 1", "Cell type 2", "Cell type 3", "Cell type 4"),
         bty = "n", fill = 1:4)
}
graphics.off()

################################

correct_idx <- c(5:8)
family_vec <- rep(c("gaussian", "poisson", "neg_binom", "curved_gaussian"), times = 3)

for(k in 1:4){
  j <- 1

  main_vec <- unlist(lapply(1:3, function(kk){
    c(paste0("Gaussian (k = ", kk, ")"),
      paste0("Poisson (k = ", kk, ")"),
      paste0("Neg. bin.\n(k = ", kk , ", est. size = ", round(res[[(k-1)*12+(kk-1)*4+3]][[j]]$fitting_param, 1), ")"),
      paste0("Curved Gaussian\n(k = ", kk , ", est. alpha = ", round(res[[(k-1)*12+(kk-1)*4+4]][[j]]$fitting_param, 1), ")"))
  }))
  main_vec[correct_idx[k]] <- paste0(main_vec[correct_idx[k]], "\n(True model)")

  png(filename = paste0("../../esvd_results/figure/experiment/Revision_writeup4_simulation_testing_", k, ".png"),
      height = 2000, width = 2750, res = 300,
      units = "px")

  par(mfrow = c(3,4), mar = c(4,4,4,0.5))
  for(i in 1:12){
    plot_prediction_against_observed(dat = res[[(k-1)*12+i]][[j]]$dat,
                                     nat_mat = res[[(k-1)*12+i]][[j]]$fit_u_mat %*% t(res[[(k-1)*12+i]][[j]]$fit_v_mat),
                                     family = family_vec[i], scalar = res[[(k-1)*12+i]][[j]]$fitting_param,
                                     missing_idx = res[[(k-1)*12+i]][[j]]$missing_idx,
                                     main = main_vec[i])
  }

  graphics.off()
}

################################
# training

correct_idx <- c(5:8)
family_vec <- rep(c("gaussian", "poisson", "neg_binom", "curved_gaussian"), times = 3)

for(k in 1:4){
  j <- 1

  main_vec <- unlist(lapply(1:3, function(kk){
    c(paste0("Gaussian (k = ", kk, ")"),
      paste0("Poisson (k = ", kk, ")"),
      paste0("Neg. bin.\n(k = ", kk , ", est. size = ", round(res[[(k-1)*12+(kk-1)*4+3]][[j]]$fitting_param, 1), ")"),
      paste0("Curved Gaussian\n(k = ", kk , ", est. alpha = ", round(res[[(k-1)*12+(kk-1)*4+4]][[j]]$fitting_param, 1), ")"))
  }))
  main_vec[correct_idx[k]] <- paste0(main_vec[correct_idx[k]], "\n(True model)")

  png(filename = paste0("../../esvd_results/figure/experiment/Revision_writeup4_simulation_training_", k, ".png"),
      height = 2000, width = 2750, res = 300,
      units = "px")

  par(mfrow = c(3,4), mar = c(4,4,4,0.5))
  for(i in 1:12){
    plot_prediction_against_observed(dat = res[[(k-1)*12+i]][[j]]$dat,
                                     nat_mat = res[[(k-1)*12+i]][[j]]$fit_u_mat %*% t(res[[(k-1)*12+i]][[j]]$fit_v_mat),
                                     family = family_vec[i], scalar = res[[(k-1)*12+i]][[j]]$fitting_param,
                                     missing_idx = c(1:prod(dim(res[[(k-1)*12+i]][[j]]$dat)))[-res[[(k-1)*12+i]][[j]]$missing_idx],
                                     main = main_vec[i])
  }

  graphics.off()
}

#################################

correct_idx <- c(5:8)

for(k in 1:4){
  main_vec <- unlist(lapply(1:3, function(kk){
    c(paste0("Gaussian (k = ", kk, ")"),
      paste0("Poisson (k = ", kk, ")"),
      paste0("Neg. bin.\n(k = ", kk , ", est. size = ", round(res[[(k-1)*12+(kk-1)*4+3]][[j]]$fitting_param, 1), ")"),
      paste0("Curved Gaussian\n(k = ", kk , ", est. alpha = ", round(res[[(k-1)*12+(kk-1)*4+4]][[j]]$fitting_param, 1), ")"))
  }))
  main_vec[correct_idx[k]] <- paste0(main_vec[correct_idx[k]], "\n(True model)")


  png(filename = paste0("../../esvd_results/figure/experiment/Revision_writeup4_simulation_embedding_", k, ".png"),
      height = 2000, width = 3500, res = 300,
      units = "px")

  par(mfrow = c(2,4), mar = c(4,4,4,0.5))

  for(i in 5:12){
    plot(res[[(k-1)*12+i]][[1]]$fit_u_mat[,1], res[[(k-1)*12+i]][[1]]$fit_u_mat[,2],
         col = rep(1:4, each = paramMat[1,"n_each"]), pch = 16,
         xlab = "Latent dimension 1", ylab = "Latent dimension 2", main = main_vec[i], asp = T)
  }

  graphics.off()
}

############################

# diagnostic

correct_idx <- c(5:8)
family_vec <- rep(c("gaussian", "poisson", "neg_binom", "curved_gaussian"), times = 3)

for(k in 1:4){
  j <- 1

  main_vec <- unlist(lapply(1:3, function(kk){
    c(paste0("Gaussian (k = ", kk, ")"),
      paste0("Poisson (k = ", kk, ")"),
      paste0("Neg. bin.\n(k = ", kk , ", est. size = ", round(res[[(k-1)*12+(kk-1)*4+3]][[j]]$fitting_param, 1), ")"),
      paste0("Curved Gaussian\n(k = ", kk , ", est. alpha = ", round(res[[(k-1)*12+(kk-1)*4+4]][[j]]$fitting_param, 1), ")"))
  }))
  main_vec[correct_idx[k]] <- paste0(main_vec[correct_idx[k]], "\n(True model)")

  png(filename = paste0("../../esvd_results/figure/experiment/Revision_writeup4_simulation_goodness_", k, ".png"),
      height = 2000, width = 2750, res = 300,
      units = "px")

  par(mfrow = c(3,4), mar = c(4,4,4,0.5))
  for(i in 1:12){
    if(i %in% c(1,5,9)){
      dat <- res[[(k-1)*12+i]][[j]]$dat
      nat_mat <- res[[(k-1)*12+i]][[j]]$fit_u_mat %*% t(res[[(k-1)*12+i]][[j]]$fit_v_mat)
      missing_idx <- res[[(k-1)*12+i]][[j]]$missing_idx
      tmp_mat <- cbind(dat[missing_idx], nat_mat[missing_idx])

      scalar <- sd(tmp_mat[,1] - tmp_mat[,2])
    } else {
      scalar <- res[[(k-1)*12+i]][[j]]$fitting_param
    }
    vec <- goodness_heuristic(dat = res[[(k-1)*12+i]][[j]]$dat,
                              nat_mat = res[[(k-1)*12+i]][[j]]$fit_u_mat %*% t(res[[(k-1)*12+i]][[j]]$fit_v_mat),
                              family = family_vec[i], scalar = scalar,
                              missing_idx = res[[(k-1)*12+i]][[j]]$missing_idx)

    plot(seq(0, 1, length.out = 21), vec, asp = T, pch = 16, xlab = "Quantile", ylab = "Empirical coverage",
         main = main_vec[i], xlim = c(0,1), ylim = c(0,1))
    lines(c(-1e2,1e2), c(-1e2,1e2), col = "red", lwd = 2, lty = 2)
  }

  graphics.off()
}

