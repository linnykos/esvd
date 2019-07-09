rm(list=ls())
load("../results/factorization_results.RData")

res_mat <- matrix(NA, 6, trials)
for(i in 1:trials){
  print(i)

  for(k in 1:6){
    cat('*')
    res_mat[k,i] <- transport::wasserstein(transport::pp(res[[1]][[i]]$dat$truth),
                                            transport::pp(res[[1]][[i]][[(k-1)*2+1]]), p = 1)
  }
}
res_mat <- res_mat[c(4,1,2,3,5,6),]

num_mat <- matrix(NA, 6, trials)
for(i in 1:trials){
  for(k in 1:6){
    num_mat[k,i] <- length(res[[1]][[i]][[2*k]]$lineages)
  }
}

# png("../figure/simulation/factorization_density.png",
#     height = 1200, width = 2000, res = 300, units = "px")
# label_vec <- c("SVD", "ICA", "t-SNE", "Our", "ZINB-WaVE", "pCMF")
# par(mfrow = c(2,3), mar = c(4, 4, 4, 0.5))
# for(i in 1:6){
#   hist(res_mat[i,], xlim = c(min(res_mat), 1), breaks = seq(min(res_mat), 1, length.out = 21), col = "gray",
#        xlab = "Kendall's tau correlation", ylab = "Count",
#        main = paste0(label_vec[i], ": (", round(median(res_mat[i,]), 2), ")"))
#   lines(rep(median(res_mat[i,]), 2), c(0, 100), col = "red", lwd = 2, lty = 2)
# }
# graphics.off()



#############################

color_func <- function(alpha = 0.2){
  c(rgb(240/255, 228/255, 66/255, alpha), #yellow
    rgb(86/255, 180/255, 233/255, alpha), #skyblue
    rgb(0/255, 158/255, 115/255, alpha), #bluish green
    rgb(0/255, 114/255, 178/255, alpha), #blue
    rgb(230/255, 159/255, 0/255, alpha), #orange
    rgb(150/255, 150/255, 150/255, alpha))
}

# start of intensive plotting function
den_list <- lapply(1:nrow(res_mat), function(i){
  density(res_mat[i,])
})

#max_val <- max(sapply(den_list, function(x){max(x$y)}))
scaling_factor <- quantile(sapply(den_list, function(x){max(x$y)}), probs = 0.3)

col_vec <- color_func(1)[c(5,2,3,1,4,6)]
text_vec <- c("eSVD", "SVD", "ICA", "t-SNE", "ZINB-WaVE", "pCMF")

png("../figure/simulation/factorization_density.png",
    height = 1800, width = 1000, res = 300, units = "px")
par(mar = c(4,0.5,4,0.5))
plot(NA, xlim = c(-.2, 1), ylim = c(0, 6.25), ylab = "",
     yaxt = "n", bty = "n", xaxt = "n", xlab = "Kendall's tau",
     main = "Estimated lineage accuracy")
axis(side = 1, at = seq(0,1,length.out = 6))
for(i in 1:nrow(res_mat)){
  lines(c(0,1), rep(nrow(res_mat) - i, 2))

  polygon(x = c(den_list[[i]]$x[1], den_list[[i]]$x, den_list[[i]]$x[length(den_list[[i]]$x)], den_list[[i]]$x[1]),
          y = (c(0, den_list[[i]]$y, 0 , 0))/scaling_factor + nrow(res_mat) - i,
          col = col_vec[i])

  med <- median(res_mat[i,])
  lines(rep(med, 2), y = c(nrow(res_mat) - i, 0), lwd = 1, lty = 2)
  points(med, y = nrow(res_mat) - i, col = "black", pch = 16, cex = 2)
  points(med, y = nrow(res_mat) - i, col = col_vec[i], pch = 16, cex = 1.5)
}
text(x = rep(0,6), y = seq(5.35,0.35,by=-1), labels = text_vec)
graphics.off()



