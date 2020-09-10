rm(list=ls())
load("../results/submission_round2b_revision/step7_additional_analyses.RData")
var_list <- ls(); var_list <- var_list[var_list != "dat_impute"]
rm(list=var_list)

####################

colorRamp_custom <- function(vec1, vec2, length){
  mat <- matrix(0, nrow = length, ncol = 3)
  for(i in 1:3){
    mat[,i] <- seq(vec1[i], vec2[i], length.out = length)
  }

  luminosity_vec <- apply(mat, 1, function(x){
    0.2126*x[1] + 0.7152*x[2] + 0.0722*x[3]
  })

  target_luminosity <- mean(c(luminosity_vec[1], luminosity_vec[length]))

  mat <- t(sapply(1:nrow(mat), function(x){
    factor <- min(c(target_luminosity/luminosity_vec[x], 1/mat[x,]))
    mat[x,] * factor
  }))

  apply(mat, 1, function(x){
    rgb(x[1], x[2], x[3])
  })
}

col_vec <- colorRamp_custom(c(0/255, 158/255, 115/255), c(213/255, 94/255, 0/255), 19)
col_vec <- c("white", col_vec)

# find the breaks
dat <- t(apply(dat_impute, 1, function(x){x/sum(x)*1000}))
dat <- dat[1:2000,]
tmp <- dat[dat != 0]
break_vec <- quantile(tmp, probs = seq(0, 1, length.out = 20))
break_vec <- c(-5, break_vec)

png("../../esvd_results/figure/experiment/heatmap.png", height = 3500, width = 2000, res = 300, units = "px")
image(.rotate(dat), breaks = break_vec, col = col_vec, asp = nrow(dat)/ncol(dat),
      axes = F)
graphics.off()

##################################



