rm(list=ls())
load("../results/factorization_esvd.RData")

correct_idx <- c(5, 18, 31, 44)
paramMat[correct_idx,]

sapply(res,length)
idx1 <- which(is.na(paramMat[,"fitting_param"]))
idx2 <- which(sapply(res, length) > 0)
idx <- intersect(idx1, idx2)
for(i in idx){print(table(sapply(res[[i]], length)))}

j <- 3
main_vec <- c("Gaussian\n(constant variance)", "Poisson",
              paste0("Negative binomial\n(true size = 50)"),
              paste0("Curved Gaussian\n(true alpha = 2)"))
family_vec <- c("gaussian", "poisson", "neg_binom", "curved_gaussian")

png(filename = "../../esvd_results/figure/experiment/Revision_writeup6_exponential_family_testing.png",
    height = 2250, width = 2250, res = 300,
    units = "px")
par(mfrow = c(2,2))
for(k in 1:4){
  i <- correct_idx[k]

  if(k %in% c(1,2)){
    nat_mat <- res[[i]][[j]]$fit$u_mat %*% t(res[[i]][[j]]$fit$v_mat)
  } else {
    nat_mat <- res[[i]][[j]]$fit[[2]]$u_mat %*% t(res[[i]][[j]]$fit[[2]]$v_mat)
  }

  if(k == 1){
    missing_idx <- res[[i]][[j]]$missing_idx
    tmp_mat <- cbind(res[[i]][[j]]$dat[missing_idx], nat_mat[missing_idx])
    scalar_val <- sd(tmp_mat[,1] - tmp_mat[,2])

  } else if(k == 3){
    scalar_val <- r_vec[2]

  } else if(k == 4){
    scalar_val <- alpha_vec[2]

  } else {
    scalar_val <- NA
  }

  plot_prediction_against_observed(dat = res[[i]][[j]]$dat,
                                   nat_mat = list(nat_mat),
                                   family = family_vec[k], scalar = scalar_val,
                                   main = main_vec[k], missing_idx_list = list(res[[i]][[j]]$missing_idx))
}

graphics.off()

####

png(filename = "../../esvd_results/figure/experiment/Revision_writeup6_exponential_family_training.png",
    height = 2250, width = 2250, res = 300,
    units = "px")
par(mfrow = c(2,2))
for(k in 1:4){
  i <- correct_idx[k]

  if(k %in% c(1,2)){
    nat_mat <- res[[i]][[j]]$fit$u_mat %*% t(res[[i]][[j]]$fit$v_mat)
  } else {
    nat_mat <- res[[i]][[j]]$fit[[2]]$u_mat %*% t(res[[i]][[j]]$fit[[2]]$v_mat)
  }

  if(k == 1){
    missing_idx <- res[[i]][[j]]$missing_idx
    tmp_mat <- cbind(res[[i]][[j]]$dat[missing_idx], nat_mat[missing_idx])
    scalar_val <- sd(tmp_mat[,1] - tmp_mat[,2])

  } else if(k == 3){
    scalar_val <- r_vec[2]

  } else if(k == 4){
    scalar_val <- alpha_vec[2]

  } else {
    scalar_val <- NA
  }

  training_idx <- c(1:prod(dim(res[[i]][[j]]$dat)))[-res[[i]][[j]]$missing_idx]

  plot_prediction_against_observed(dat = res[[i]][[j]]$dat,
                                   nat_mat = list(nat_mat),
                                   family = family_vec[k], scalar = scalar_val,
                                   main = main_vec[k], missing_idx_list = list(training_idx))
}

graphics.off()

######################

j <- 3
png(filename = "../../esvd_results/figure/experiment/Revision_writeup6_exponential_family_fit.png",
    height = 2250, width = 2250, res = 300,
    units = "px")
pos_vec <- c("topleft", "topleft", "topright", "topright")
par(mfrow = c(2,2))
for(k in 1:4){
  i <- correct_idx[k]

  if(k %in% c(1,2)){
    u_mat <- res[[i]][[j]]$fit$u_mat
  } else {
    u_mat <- res[[i]][[j]]$fit[[2]]$u_mat
  }

  i <- correct_idx[k]
  plot(u_mat[,1], u_mat[,2], col = rep(1:4, each = paramMat[1,"n_each"]), pch = 16,
       xlab = "Latent dimension 1", ylab = "Latent dimension 2", main =
         paste0("Estimated cell latent positions\n(", main_vec[k], ")"), asp = T)
  legend(pos_vec[k], c("Cell type 1", "Cell type 2", "Cell type 3", "Cell type 4"),
         bty = "n", fill = 1:4)
}
graphics.off()
