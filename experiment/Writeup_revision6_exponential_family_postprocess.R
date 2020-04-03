rm(list=ls())
load("../results/factorization_exponential_families.RData")

correct_idx <- c(5, 18, 31, 44)
paramMat[correct_idx,]

sapply(res,length)
idx1 <- which(is.na(paramMat[,"fitting_param"]))
idx2 <- which(sapply(res, length) > 0)
nusiance_idx <- intersect(idx1, idx2)
for(i in idx){print(table(sapply(res[[i]], length)))}

# some basic cleanup
for(i in 1:length(res)){
  len_vec <- sapply(res[[i]], length)
  res[[i]] <- res[[i]][which(len_vec > 1)]
}

# select which nuisance parameter is best, among all the neg binom and curved gaussians
selected_idx <- vector("list", length(res))

for(i in nusiance_idx){
  print(i)

  # determine among this parameter setting, how each of 100 trials behave
  tmp_vec <- sapply(res[[i]], function(x){

    # if i is odd, then it's a negative binomial. if i is even, then it's a curved gaussian
    if(i %% 2 == 1) {
      scalar_vec <- r_vec; family_val <- "neg_binom"
    } else {
      scalar_vec <- alpha_vec; family_val <- "curved_gaussian"
    }

    quality_vec <- sapply(1:length(x$fit), function(j){
      nat_mat <- x$fit[[j]]$u_mat %*% t(x$fit[[j]]$v_mat)
      eSVD::plot_prediction_against_observed(dat = x$dat, nat_mat_list = list(nat_mat),
                                       family = family_val, missing_idx_list = list(x$missing_idx),
                                       scalar = scalar_vec[j], plot = F)
    })


    which.min(abs(quality_vec - 45))
  })

  # tmp_vec <- factor(as.numeric(tmp_vec), levels = c(1,2,3))
  # selected_idx[[i]] <- table(tmp_vec)
  selected_idx[[i]] <- tmp_vec
}

#########################

j <- 1
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

j <- 1
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

###############################

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
