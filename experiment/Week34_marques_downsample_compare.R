rm(list=ls())
load("../experiment/Week33_downsample.RData")

downsample_viper_vec <- numeric(9)
dropout_viper_vec <- numeric(9)
make_plot <- TRUE

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

