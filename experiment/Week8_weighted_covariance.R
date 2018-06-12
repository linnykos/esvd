rm(list=ls())

# set up parameters
adj <- matrix(c(0, 0.1, -.05,
                0.15,-.1,-.1,
                -.25,-.15,0.05), 3, 3, byrow = T)
tmp <- svd(adj)
u_center <- t(tmp$u %*% diag(sqrt(tmp$d)))
v_center <- t(tmp$v %*% diag(sqrt(tmp$d)))

t(u_center) %*% v_center

u_num <- c(20, 10, 80)
u_label <- unlist(lapply(1:3, function(x){rep(x, u_num[x])}))
v_num <- c(30, 40, 100)
v_label <- unlist(lapply(1:3, function(x){rep(x, v_num[x])}))

# generate matrices
set.seed(10)
u_sig <- 0.03
u_dat <- do.call(rbind, lapply(1:3, function(x){
  MASS::mvrnorm(n = u_num[x], mu = u_center[,x], Sigma = u_sig*diag(3))
}))
v_sig <- 0.03
v_dat <- do.call(rbind, lapply(1:3, function(x){
  MASS::mvrnorm(n = v_num[x], mu = v_center[,x], Sigma = v_sig*diag(3))
}))

mean_dat <- u_dat %*% t(v_dat)
n <- nrow(mean_dat)
d <- ncol(mean_dat)
dat <- matrix(0, n, d)
setting <- matrix(0, n, d)
nodropout_dat <- matrix(0, n, d)

dropout_create <- function(a, b){
  function(x){ 1/(1+exp(-(a+b*x))) }
}
dropout_coef <- c(-1, 5)
dropout_function <- dropout_create(dropout_coef[1], dropout_coef[2])

set.seed(10)
for(i in 1:n){
  for(j in 1:d){
    if(mean_dat[i,j] <= 0) {
      dat[i,j] <- 0
      setting[i,j] <- 0
      nodropout_dat[i,j] <- 0
    } else {
      val <- rexp(1, rate = 1/mean_dat[i,j])
      bool <- rbinom(1, 1, prob = dropout_function(val))
      setting[i,j] <- bool + 1
      dat[i,j] <- bool*val
      nodropout_dat[i,j] <- mean_dat[i,j]
    }
  }
}

#############

weighted_correlation_paired <- function(vec1, vec2, dropout_coef, kappa = 0.95){
  stopifnot(length(vec1) == length(vec2))
  l <- length(vec1)
  dropout_function <- dropout_create(dropout_coef[1], dropout_coef[2])

  prob1 <- 1-sapply(vec1, dropout_function) #probability of observing dropout
  prob2 <- 1-sapply(vec2, dropout_function)

  weight_vec <- kappa*sqrt((1-prob1)*(1-prob2)) + (1-kappa)
  boot::corr(cbind(vec1, vec2), w = weight_vec)
}

weighted_correlation <- function(dat, dropout_coef, kappa = 0.95){
  l <- ncol(dat)
  mat <- matrix(0, l, l)

  for(i in 1:l){
    for(j in i:l){
      mat[i,j] <- weighted_correlation_paired(dat[,i], dat[,j], dropout_coef, kappa)
      mat[j,i] <- mat[i,j]
    }
  }

  mat
}

weighted_gene_cor <- weighted_correlation(dat, dropout_coef, 0.95)
gene_cor <- cor(dat)
true_cor <- cor(nodropout_dat)

col_idx <- sapply(1:(max(v_label)-1), function(x){
  length(which(v_label <= x))/length(v_label)
})

col_vec2 <- colorRampPalette(c(rgb(0.584, 0.858, 0.564), rgb(0.803, 0.156, 0.211)))(19)
break_vec <- quantile(c(as.numeric(true_cor), as.numeric(weighted_gene_cor),
                        as.numeric(gene_cor)), probs = seq(0, 1, length.out = 20))
.rotate <- function(mat){t(mat)[,nrow(mat):1]}

png(paste0("../figure/experiment/8_simulated_weighted_gene_correlation.png"), height = 2400, width = 7200, res = 300, units = "px")
par(mfrow = c(1,3))
image(.rotate(true_cor), breaks = break_vec, col = col_vec2, asp = T, axes = F,
      main = "True correlation", cex.main = 3)

for(j in col_idx){
  lines(c(0,1), rep(1-j, 2), lwd = 6, col = "white")
  lines(c(0,1), rep(1-j, 2), lwd = 6, lty = 2)
}
for(j in col_idx){
  lines(rep(j, 2), c(0,1), lwd = 6, col = "white")
  lines(rep(j, 2), c(0,1), lwd = 6, lty = 2)
}

image(.rotate(gene_cor), breaks = break_vec, col = col_vec2, asp = T, axes = F,
      main = "Estimated correlation", cex.main = 3)

for(j in col_idx){
  lines(c(0,1), rep(1-j, 2), lwd = 6, col = "white")
  lines(c(0,1), rep(1-j, 2), lwd = 6, lty = 2)
}
for(j in col_idx){
  lines(rep(j, 2), c(0,1), lwd = 6, col = "white")
  lines(rep(j, 2), c(0,1), lwd = 6, lty = 2)
}

image(.rotate(weighted_gene_cor), breaks = break_vec, col = col_vec2, asp = T, axes = F,
      main = "Weighted estimated correlation", cex.main = 3)

for(j in col_idx){
  lines(c(0,1), rep(1-j, 2), lwd = 6, col = "white")
  lines(c(0,1), rep(1-j, 2), lwd = 6, lty = 2)
}
for(j in col_idx){
  lines(rep(j, 2), c(0,1), lwd = 6, col = "white")
  lines(rep(j, 2), c(0,1), lwd = 6, lty = 2)
}

graphics.off()

############

## see a plot of comparison

png(paste0("../figure/experiment/8_simulated_weighted_correlation_comparison.png"), height = 1350, width = 2400, res = 300, units = "px")
par(mfrow = c(1,2))

#plot(as.numeric(true_cor), as.numeric(gene_cor), asp = T, pch = 16, cex.main = 3)
#lines(c(-2,2), c(-2,2), col = "red", lwd = 2)
plot(as.numeric(true_cor), as.numeric(weighted_gene_cor), asp = T, pch = 16, cex.main = 3,
     col = rgb(0,0,0,0.02), xlab = "True correlation value",
     ylab = "Weighted correlation value")
lines(c(-2,2), c(-2,2), col = "red", lwd = 2)
plot(as.numeric(gene_cor), as.numeric(weighted_gene_cor), asp = T, pch = 16, cex.main = 3,
     col = rgb(0,0,0,0.02), xlab = "Pearson correlation value",
     ylab = "Weighted correlation value")
lines(c(-2,2), c(-2,2), col = "red", lwd = 2)

graphics.off()

### # find the instance where weighted correlation makes the most difference

color_correlation_func <- function(vec1, vec2, dropout_coef, kappa = 0.95){
  dropout_function <- dropout_create(dropout_coef[1], dropout_coef[2])

  prob1 <- 1-sapply(vec1, dropout_function) #probability of observing dropout
  prob2 <- 1-sapply(vec2, dropout_function)

  kappa*(1-prob1)*(1-prob2) + (1-kappa)
}

break_vec <- seq(0, 1, length.out = 19)

png(paste0("../figure/experiment/8_simulated_gene_scatterplot.png"), height = 900, width = 2400, res = 300, units = "px")
par(mfrow = c(1,3))

# largest difference
val <- max(abs(gene_cor - weighted_gene_cor))
idx <- which(abs(gene_cor - weighted_gene_cor) == val, arr.ind = T)
vec <- color_correlation_func(dat[,idx[1]], dat[,idx[2]], c(-2,30))
col_vec <- sapply(vec, function(x){col_vec2[which.min(abs(break_vec - x))]})
plot(dat[,idx[1]], dat[,idx[2]], pch = 16, col = col_vec, cex = 1.5,
     xlab = paste0("Gene ", idx[1]), ylab = paste0("Gene ", idx[2]),
     main = paste0("Largest difference,\nPearson: ", round(cor(dat[,idx[1]], dat[,idx[2]]), 2)))

# largest difference for negative pearson
tmp <- which(gene_cor < 0)
val <- max(abs(gene_cor[tmp] - weighted_gene_cor[tmp]))
idx <- which(abs(gene_cor - weighted_gene_cor) == val, arr.ind = T)
vec <- color_correlation_func(dat[,idx[1]], dat[,idx[2]], c(-2,30))
col_vec <- sapply(vec, function(x){col_vec2[which.min(abs(break_vec - x))]})
plot(dat[,idx[1]], dat[,idx[2]], pch = 16, col = col_vec, cex = 1.5,
     xlab = paste0("Gene ", idx[1]), ylab = paste0("Gene ", idx[2]),
     main = paste0("Largest difference for negatives,\nPearson: ", round(cor(dat[,idx[1]], dat[,idx[2]]), 2)))

# smallest difference
tmp <- which(gene_cor != 1)
val <- min(abs(gene_cor[tmp] - weighted_gene_cor[tmp]))
idx <- which(abs(gene_cor - weighted_gene_cor) == val, arr.ind = T)
vec <- color_correlation_func(dat[,idx[1]], dat[,idx[2]], c(-2,30))
col_vec <- sapply(vec, function(x){col_vec2[which.min(abs(break_vec - x))]})
plot(dat[,idx[1]], dat[,idx[2]], pch = 16, col = col_vec, cex = 1.5,
     xlab = paste0("Gene ", idx[1]), ylab = paste0("Gene ", idx[2]),
     main = paste0("Smallest difference,\nPearson: ", round(cor(dat[,idx[1]], dat[,idx[2]]), 2)))

graphics.off()

##########################

# dropout curves
set.seed(10)
fit_coef <- estimate_dropout(dat, threshold_quant_degree = 1, threshold_quant_logistic = 0, k = 3)
fit_coef
est_dropout_func <- dropout_create(fit_coef[1], fit_coef[2])

png(paste0("../figure/experiment/8_dropout_function.png"), height = 1200, width = 2400, res = 300, units = "px")
par(mfrow = c(1,2))

x_vec <- seq(0, max(dat), length.out = 100)
val1 <- sapply(x_vec, dropout_function)
val2 <- sapply(x_vec, est_dropout_func)
plot(x_vec, val1, pch = 16, xlab = "Observed value", ylab = "Probability of no dropout")
points(x_vec, val2, col = rgb(0.803, 0.156, 0.211), pch = 16)

legend("bottomright", c("True", "Estimated"),
       bty="n", fill=c("black", rgb(0.803, 0.156, 0.211)))

quant_vec <- seq(0, 1, length.out = 100)
idx <- which(dat != 0)
x_vec <- sapply(quant_vec, function(x){quantile(dat[idx], probs = x)})
val1 <- sapply(x_vec, dropout_function)
val2 <- sapply(x_vec, est_dropout_func)
plot(quant_vec, val1, pch = 16, xlab = "Quantile of observed value", ylab = "Probability of no dropout")
points(quant_vec, val2, col = rgb(0.803, 0.156, 0.211), pch = 16)

legend("topleft", c("True", "Estimated"),
       bty="n", fill=c("black", rgb(0.803, 0.156, 0.211)))

graphics.off()

##############################

# see if number of zeros is correlated with PCA
# let's do genes

