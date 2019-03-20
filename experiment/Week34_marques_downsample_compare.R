rm(list=ls())
load("../experiment/Week33_downsample.RData")

downsample_viper_vec <- numeric(9)
dropout_viper_vec <- numeric(9)
make_plot <- FALSE

for(i in 1:9){
  tmp <- t(viper_list[[i]]$imputed)

  downsample_idx <- downsample_list[[i]]$downsample_idx
  downsample_viper_vec[i] <- sum(abs(tmp[downsample_idx] - dat[downsample_idx]))/length(downsample_idx)

  if(make_plot){
    plot(dat[downsample_idx], tmp[downsample_idx], asp = T, pch = 16,
         col = rgb(0,0,0,0.2), main = "VIPER Downsample")
    lines(c(0, 10000), c(0, 10000), col = "red", lwd = 2)
  }

  dropout_idx <- downsample_list[[i]]$dropout_idx
  dropout_viper_vec[i] <- cor(tmp[dropout_idx], dat[dropout_idx])

  if(make_plot){
    plot(dat[dropout_idx], tmp[dropout_idx], asp = T, pch = 16,
         col = rgb(0,0,0,0.1), main = "VIPER Dropout")
    lines(c(0, 10000), c(0, 10000), col = "red", lwd = 2)
  }
}

#############################

downsample_scimpute_vec <- numeric(9)
dropout_scimpute_vec <- numeric(9)
make_plot <- TRUE

for(i in 1:9){
  tmp <- t(read.csv(paste0("../experiment/Week33_scimpute_", i, "scimpute_count.csv")))[-1,]

  downsample_idx <- downsample_list[[i]]$downsample_idx
  downsample_scimpute_vec[i] <- sum(abs(tmp[downsample_idx] - dat[downsample_idx]))/length(downsample_idx)

  if(make_plot){
    plot(dat[downsample_idx], tmp[downsample_idx], asp = T, pch = 16,
         col = rgb(0,0,0,0.2), main = "scImpute Downsample")
    lines(c(0, 10000), c(0, 10000), col = "red", lwd = 2)
  }

  dropout_idx <- downsample_list[[i]]$dropout_idx
  dropout_scimpute_vec[i] <- cor(tmp[dropout_idx], dat[dropout_idx])

  if(make_plot){
    plot(dat[dropout_idx], tmp[dropout_idx], asp = T, pch = 16,
         col = rgb(0,0,0,0.1), main = "scImpute Dropout")
    lines(c(0, 10000), c(0, 10000), col = "red", lwd = 2)
  }
}

x_lim <- range(c(dropout_viper_vec, dropout_scimpute_vec)); y_lim <- range(c(downsample_scimpute_vec, downsample_viper_vec))

png("../figure/experiment/Week34_marques_downsampling.png", height = 1200, width = 1200, res = 300, units = "px")
plot(NA, xlim = x_lim, ylim = y_lim, main = "Downsample results", xlab = "Correlation (drop-out entries)",
     ylab = "L1 loss (down-sampling entries)")

points(dropout_viper_vec, downsample_viper_vec, pch = 16, cex = 2)
points(dropout_scimpute_vec, downsample_scimpute_vec, pch = 16, cex = 2, col = "red")


legend("topleft",  c("VIPER", "scImpute"),
       fill = c(1,2), bty = "n", cex = 0.8)

graphics.off()

##########

i <- 5
tmp <- t(viper_list[[i]]$imputed)
downsample_idx <- downsample_list[[i]]$downsample_idx
downsample_viper_vec[i] <- sum(abs(tmp[downsample_idx] - dat[downsample_idx]))/length(downsample_idx)
dropout_idx <- downsample_list[[i]]$dropout_idx
dropout_viper_vec[i] <- cor(tmp[dropout_idx], dat[dropout_idx])

png("../figure/experiment/Week34_marques_downsampling_VIPER.png", height = 1200, width = 2000, res = 300, units = "px")
par(mfrow = c(1,2))
plot(dat[downsample_idx], tmp[downsample_idx], asp = T, pch = 16,
     col = rgb(0,0,0,0.2), main = paste0("VIPER Downsample\nL1 = ", round(downsample_viper_vec[i], 2)), xlim = c(0,60), ylim = c(0,60),
     xlab = "Masked value (observed)", ylab = "Imputed value (estimated)")
lines(c(0, 10000), c(0, 10000), col = "red", lwd = 2)

plot(dat[dropout_idx], tmp[dropout_idx], asp = T, pch = 16,
     col = rgb(0,0,0,0.1), main = paste0("VIPER Dropout\nCorrelation = ", round(dropout_viper_vec[i], 2)), xlim = c(0, 200), ylim = c(0, 200),
     xlab = "Masked value (observed)", ylab = "Imputed value (estimated)")
lines(c(0, 10000), c(0, 10000), col = "red", lwd = 2)
graphics.off()

##

tmp <- t(read.csv(paste0("../experiment/Week33_scimpute_", i, "scimpute_count.csv")))[-1,]

downsample_idx <- downsample_list[[i]]$downsample_idx
downsample_scimpute_vec[i] <- sum(abs(tmp[downsample_idx] - dat[downsample_idx]))/length(downsample_idx)
dropout_idx <- downsample_list[[i]]$dropout_idx
dropout_scimpute_vec[i] <- cor(tmp[dropout_idx], dat[dropout_idx])

png("../figure/experiment/Week34_marques_downsampling_scImpute.png", height = 1200, width = 2000, res = 300, units = "px")
par(mfrow = c(1,2))
plot(dat[downsample_idx], tmp[downsample_idx], asp = T, pch = 16,
     col = rgb(0,0,0,0.2), main = paste0("scImpute Downsample\nL1 = ", round(downsample_scimpute_vec[i], 2)), xlim = c(0,60), ylim = c(0,60),
     xlab = "Masked value (observed)", ylab = "Imputed value (estimated)")
lines(c(0, 10000), c(0, 10000), col = "red", lwd = 2)

plot(dat[dropout_idx], tmp[dropout_idx], asp = T, pch = 16,
     col = rgb(0,0,0,0.1), main = paste0("scImpute Dropout\nCorrelation = ", round(dropout_scimpute_vec[i], 2)), xlim = c(0, 200), ylim = c(0, 200),
     xlab = "Masked value (observed)", ylab = "Imputed value (estimated)")
lines(c(0, 10000), c(0, 10000), col = "red", lwd = 2)
graphics.off()

