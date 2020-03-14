rm(list=ls())
load("../results/lingxue_analysis.RData")

#fit_all_list <- fit_all_list[1:6]
lapply(fit_all_list, function(x){x$neg_bin_param})
lapply(fit_all_list, function(x){x$curved_gaussian_param})


# 1,2,3,4,5
# which dataset to investigate?
for(dat_num in 1:length(fit_all_list)){
  label_vec <- preprocessing_list[[dat_num]]$label_vec
  label_num <- as.numeric(label_vec)

  main_vec <- c("Gaussian", "Poisson", "Exponential", "Negative binomial", "Curved Gaussian")
  main_vec[4] <- paste0(main_vec[4], ": Size=", round(fit_all_list[[dat_num]]$neg_bin_param[length(fit_all_list[[dat_num]]$neg_bin_param)], 2))
  main_vec[5] <- paste0(main_vec[5], ": Alpha=", round(fit_all_list[[dat_num]]$curved_gaussian_param[length(fit_all_list[[dat_num]]$curved_gaussian_param)],2))
  idx_vec <- seq(1, 9, by = 2)

  png(filename = paste0("../figure/experiment/Revision_writeup4_lingxue_embedding_", dat_num, "_k", k, ".png"),
      height = 2000, width = 2500, res = 300,
      units = "px")
  par(mfrow = c(2:3))
  for(i in 1:5){
    plot(fit_all_list[[dat_num]][[idx_vec[i]]]$u_mat[,1], fit_all_list[[dat_num]][[idx_vec[i]]]$u_mat[,2],
         asp = T, col = c(1:5)[label_num],
         xlab = "Latent dimension 1", ylab = "Latent dimension 2",
         pch = 16, main = paste0("Dataset ", dat_num, main_vec[i]))
  }
  graphics.off()
}

###############################

for(dat_num in 1:length(fit_all_list)){
  dat_impute <- preprocessing_list[[dat_num]]$dat_impute
  set.seed(10)
  missing_idx <- construct_missing_values(n = nrow(dat_impute), p = ncol(dat_impute), num_val = 2)

  family_vec <- c("gaussian", "poisson", "exponential", "neg_binom", "curved_gaussian")
  main_vec <- c("Gaussian", "Poisson", "Exponential", "Negative binomial", "Curved Gaussian")
  main_vec[4] <- paste0(main_vec[4], ": Size=", round(fit_all_list[[dat_num]]$neg_bin_param[length(fit_all_list[[dat_num]]$neg_bin_param)], 2))
  main_vec[5] <- paste0(main_vec[5], ": Alpha=", round(fit_all_list[[dat_num]]$curved_gaussian_param[length(fit_all_list[[dat_num]]$curved_gaussian_param)],2))
  idx_vec <- seq(1, 9, by = 2)

  png(filename = paste0("../figure/experiment/Revision_writeup4_lingxue_missing_", dat_num, "_k", k, ".png"),
      height = 2000, width = 2500, res = 300,
      units = "px")
  par(mfrow = c(2:3))
  for(i in 1:5){
    nat_mat <- fit_all_list[[dat_num]][[idx_vec[i]]]$u_mat %*% t(fit_all_list[[dat_num]][[idx_vec[i]]]$v_mat)
    scalar <- ifelse(family_vec[i] == "curved_gaussian", fit_all_list[[dat_num]]$curved_gaussian_param[length(fit_all_list[[dat_num]]$curved_gaussian_param)],
                     fit_all_list[[dat_num]]$neg_bin_param[length(fit_all_list[[dat_num]]$neg_bin_param)])
    pred_mat <- compute_mean(nat_mat, family = family_vec[i], scalar = scalar)

    plot_mat <- cbind(dat_impute[missing_idx], pred_mat[missing_idx])

    seq_val <- seq(0, 4000, length.out = 500)
    if(i == 1){
      sd_val <- stats::sd(plot_mat[,1] - plot_mat[,2])

      y_bot <- sapply(seq_val, function(x){
        stats::qnorm(0.1, mean = x, sd = sd_val)
      })
      y_top <- sapply(seq_val, function(x){
        stats::qnorm(0.9, mean = x, sd = sd_val)
      })
    } else if(i == 2){
      y_bot <- sapply(seq_val, function(x){
        stats::qpois(0.1, lambda = x)
      })
      y_top <- sapply(seq_val, function(x){
        stats::qpois(0.9, lambda = x)
      })
    } else if(i == 3){
      y_bot <- sapply(seq_val, function(x){
        stats::qexp(0.1, rate = 1/x)
      })
      y_top <- sapply(seq_val, function(x){
        stats::qexp(0.9, rate = 1/x)
      })
    } else if(i == 4){
      y_bot <- sapply(seq_val, function(x){
        stats::qnbinom(0.1, size = scalar, prob = scalar/(scalar+x))
      })
      y_top <- sapply(seq_val, function(x){
        stats::qnbinom(0.9, size = scalar, prob = scalar/(scalar+x))
      })
    } else {
      y_bot <- sapply(seq_val, function(x){
        stats::qnorm(0.1, mean = x, sd = x/scalar)
      })
      y_top <- sapply(seq_val, function(x){
        stats::qnorm(0.9, mean = x, sd = x/scalar)
      })
    }

    plot(NA, asp = T, main = paste0("Dataset ", dat_num, "\n", main_vec[i]), pch = 16,
         xlim = range(plot_mat), ylim = range(plot_mat),
         xlab = "Predicted missing value", ylab = "Observed value")

    polygon(c(seq_val, rev(seq_val)), c(y_top, rev(y_bot)), col = rgb(1,0,0,0.2),
            border = NA, density = 30, angle = -45)
    points(plot_mat[,2], plot_mat[,1], pch = 16, col = rgb(0,0,0,0.2))

    lines(rep(0,2), c(-1e10,1e10), col = "red", lwd = 1)
    lines(c(-1e10,1e10), rep(0,2), col = "red", lwd = 1)
    lines(c(-1e10,1e10), c(-1e10,1e10), col = "red", lwd = 2)

    lines(seq_val, y_bot, col = "red", lty = 2, lwd = 2)
    lines(seq_val, y_top, col = "red", lty = 2, lwd = 2)

    pca_res <- stats::prcomp(plot_mat, center = F, scale = F)
    lines(c(0, 1e6), c(0, 1e6*pca_res$rotation[2,1]/pca_res$rotation[1,1]), col = "blue", lwd = 2, lty = 2)
  }

  graphics.off()
}

for(dat_num in 1:length(fit_all_list)){
  dat_impute <- preprocessing_list[[dat_num]]$dat_impute
  set.seed(10)
  missing_idx <- construct_missing_values(n = nrow(dat_impute), p = ncol(dat_impute), num_val = 2)

  family_vec <- c("gaussian", "poisson", "exponential", "neg_binom", "curved_gaussian")
  main_vec <- c("Gaussian", "Poisson", "Exponential", "Negative binomial", "Curved Gaussian")
  main_vec[4] <- paste0(main_vec[4], ": Size=", round(fit_all_list[[dat_num]]$neg_bin_param[length(fit_all_list[[dat_num]]$neg_bin_param)], 2))
  main_vec[5] <- paste0(main_vec[5], ": Alpha=", round(fit_all_list[[dat_num]]$curved_gaussian_param[length(fit_all_list[[dat_num]]$curved_gaussian_param)],2))
  idx_vec <- seq(2, 10, by = 2)

  png(filename = paste0("../figure/experiment/Revision_writeup4_lingxue_training_", dat_num, "_k", k, ".png"),
      height = 2000, width = 2500, res = 300,
      units = "px")
  par(mfrow = c(2:3))
  for(i in 1:5){
    nat_mat <- fit_all_list[[dat_num]][[idx_vec[i]]]$u_mat %*% t(fit_all_list[[dat_num]][[idx_vec[i]]]$v_mat)
    scalar <- ifelse(family_vec[i] == "curved_gaussian", fit_all_list[[dat_num]]$curved_gaussian_param[length(fit_all_list[[dat_num]]$curved_gaussian_param)],
                     fit_all_list[[dat_num]]$neg_bin_param[length(fit_all_list[[dat_num]]$neg_bin_param)])
    pred_mat <- compute_mean(nat_mat, family = family_vec[i], scalar = scalar)

    plot_mat <- cbind(dat_impute[missing_idx], pred_mat[missing_idx])

    seq_val <- seq(0, 4000, length.out = 500)
    if(i == 1){
      sd_val <- stats::sd(plot_mat[,1] - plot_mat[,2])

      y_bot <- sapply(seq_val, function(x){
        stats::qnorm(0.1, mean = x, sd = sd_val)
      })
      y_top <- sapply(seq_val, function(x){
        stats::qnorm(0.9, mean = x, sd = sd_val)
      })
    } else if(i == 2){
      y_bot <- sapply(seq_val, function(x){
        stats::qpois(0.1, lambda = x)
      })
      y_top <- sapply(seq_val, function(x){
        stats::qpois(0.9, lambda = x)
      })
    } else if(i == 3){
      y_bot <- sapply(seq_val, function(x){
        stats::qexp(0.1, rate = 1/x)
      })
      y_top <- sapply(seq_val, function(x){
        stats::qexp(0.9, rate = 1/x)
      })
    } else if(i == 4){
      y_bot <- sapply(seq_val, function(x){
        stats::qnbinom(0.1, size = scalar, prob = scalar/(scalar+x))
      })
      y_top <- sapply(seq_val, function(x){
        stats::qnbinom(0.9, size = scalar, prob = scalar/(scalar+x))
      })
    } else {
      y_bot <- sapply(seq_val, function(x){
        stats::qnorm(0.1, mean = x, sd = x/scalar)
      })
      y_top <- sapply(seq_val, function(x){
        stats::qnorm(0.9, mean = x, sd = x/scalar)
      })
    }

    plot(NA, asp = T, main = paste0("Dataset ", dat_num, "\n", main_vec[i]), pch = 16,
         xlim = range(plot_mat), ylim = range(plot_mat),
         xlab = "Predicted missing value", ylab = "Observed value")

    polygon(c(seq_val, rev(seq_val)), c(y_top, rev(y_bot)), col = rgb(1,0,0,0.2),
            border = NA, density = 30, angle = -45)
    points(plot_mat[,2], plot_mat[,1], pch = 16, col = rgb(0,0,0,0.2))

    lines(rep(0,2), c(-1e10,1e10), col = "red", lwd = 1)
    lines(c(-1e10,1e10), rep(0,2), col = "red", lwd = 1)
    lines(c(-1e10,1e10), c(-1e10,1e10), col = "red", lwd = 2)

    lines(seq_val, y_bot, col = "red", lty = 2, lwd = 2)
    lines(seq_val, y_top, col = "red", lty = 2, lwd = 2)

    pca_res <- stats::prcomp(plot_mat, center = F, scale = F)
    lines(c(0, 1e6), c(0, 1e6*pca_res$rotation[2,1]/pca_res$rotation[1,1]), col = "blue", lwd = 2, lty = 2)
  }

  graphics.off()
}


