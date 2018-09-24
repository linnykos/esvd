set.seed(10)
.rotate <- function(mat){t(mat)[,nrow(mat):1]}

vec <- c(stats::rexp(n = 2000),
         stats::rgamma(n = 2000, shape = 2, scale = 2) + 3,
         stats::rnorm(n = 1995, mean = 17))

hist(vec, breaks = 100, col = "gray")

vec <- sort(vec, decreasing = T)

# make a 5-block design
design <- matrix(0, d, d)
for(k in 1:5){
  for(i in 1:(5-k+1)){
    design[((i-1)*(d/5)+1):(i*d/5), (((i-1)*(d/5)+1):(i*d/5)) + (k-1)*(d/5)] <- 10*i + (5-k)*100
    design[(((i-1)*(d/5)+1):(i*d/5)) + (k-1)*(d/5), ((i-1)*(d/5)+1):(i*d/5)] <- 10*i + (5-k)*100
  }
}

for(i in 1:d){
  for(j in 1:i){
    design[i,j] <- stats::rnorm(1, mean = design[i,j], sd = 0.1)
  }
}
design_vec <- design[lower.tri(design, diag = F)]
idx_order <- order(design_vec, decreasing = T)

d <- 110
mat <- matrix(0, d, d)
mat_vec <- rep(NA, n)
for(i in 1:n){
  mat_vec[idx_order[i]] <- vec[i]
}
mat[lower.tri(mat, diag = F)] <- mat_vec
mat <- mat + t(mat)

image(mat)

eigen_res <- eigen(mat)
plot(eigen_res$values)

# mat is the "true" matrix of means
obs <- rep(NA, n)
for(i in 1:n){
  obs[i] <- stats::rnorm(1, mean = mat_vec[i], sd = 0.1)
}
hist(obs, breaks = 100, col = "gray")

png("../figure/theory/P2_matrix.png", height = 1000, width = 2400, res = 300, units = "px")
par(mfrow = c(1,2))
par(mar = rep(0.5,4))
image(.rotate(mat), asp = T, xlab = "", ylab = "", axes = F)

par(mar = c(4,4,4,0))
hist(obs, breaks = 100, col = "gray", xlab = "Observed values", ylab = "Count",
     main = "Histogram of observed values")
graphics.off()

png("../figure/theory/P2_eigen.png", height = 1200, width = 2400, res = 300, units = "px")
par(mfrow = c(1,2))
plot(eigen_res$values, pch = 16, xlab = "Index",
     ylab = "Value", main = "Population eigenvalue")
lines(c(-1e6, 1e6), rep(0, 2), col = "red", lwd = 2)

ylim = range(eigen_res$vectors)
plot(NA, xlim = c(1, nrow(eigen_res$vectors)), ylim = ylim,
     main = "Top 5 (positive)\npopulation eigenvectors", xlab = "Index", ylab = "Value")
idx <- 5
points(eigen_res$vectors[,idx], pch = 16, col = 5)
lines(eigen_res$vectors[,idx], lwd = 2, col = 5)

for(i in 1:4){
  points(eigen_res$vectors[,i], pch = 16, col = i)
  lines(eigen_res$vectors[,i], lwd = 2, col = i)
}

graphics.off()

#########################################

.decompose <- function(mat, k = 1){
  eigen_res <- eigen(mat)
  idx <- order(abs(eigen_res$values), decreasing = T)

  if(k == 1) {
    diag_mat <- as.matrix(eigen_res$values[idx[1:k]])
  } else {
    diag_mat <- diag(eigen_res$values[idx[1:k]])
  }

  recon <- eigen_res$vectors[,idx[1:k],drop = F] %*% diag_mat %*% t(eigen_res$vectors[,idx[1:k], drop = F])
  recon_vec <- recon[lower.tri(recon, diag = F)]

  vec <- rep(NA, length(recon_vec))
  for(i in 1:length(recon_vec)){
    vec[i] <- stats::rnorm(1, mean = recon_vec[i], sd = 0.1)
  }

  vec
}

seq_breaks <- seq(-5, 30, length.out = 100)
k_vec <- c(1:4, 5, 10, 20, 110)

tmp_list <- lapply(k_vec, function(k){
  set.seed(10)
  .decompose(mat, k = k)
})

png("../figure/theory/P2_hist.png", height = 2500, width = 2000, res = 300, units = "px")
par(mfrow = c(4,2), mar = c(4,4,4,0.5))
for(i in 1:8){
  hist(tmp_list[[i]], breaks = seq_breaks, col = "gray",
       main = paste0("Using rank ", k_vec[i], " approximation"),
       xlab = "Value", ylab = "Count")
}
graphics.off()



