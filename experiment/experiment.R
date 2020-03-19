rm(list=ls())

load("../results/factorization_esvd.RData")
seq_length <- 21

i <- 19
j <- 1
dat <- res[[i]][[j]]$dat
nat_mat <- res[[i]][[j]]$fit_u_mat %*% t(res[[i]][[j]]$fit_v_mat)
family <- "neg_binom"
missing_idx <- res[[i]][[j]]$missing_idx

pred_mat <- compute_mean(nat_mat, family = family, scalar = res[[i]][[j]]$fitting_param)

tmp_mat <- cbind(dat[missing_idx], pred_mat[missing_idx])

width_seq <- seq(0, 1, length = seq_length)

vec <- sapply(width_seq, function(width){
  interval_mat <- sapply(tmp_mat[,2], function(x){
    .compute_prediction_interval_from_mean(x, family = family, width = width, scalar = res[[i]][[j]]$fitting_param)
  })

  sum(sapply(1:nrow(tmp_mat), function(x){
    interval_mat[1,x] <= tmp_mat[x,1] & interval_mat[2,x] >= tmp_mat[x,1]
  }))/nrow(tmp_mat)
})

plot(width_seq, vec, asp = T, xlim = c(0,1), ylim = c(0,1), pch = 16)
lines(c(-1e3,1e3), c(-1e3,1e3), col = "red", lty = 2, lwd = 2)

plot_prediction_against_observed(dat, nat_mat, family = "neg_binom", missing_idx = missing_idx, main = "Test",
                                 scalar = res[[i]][[j]]$fitting_param)


