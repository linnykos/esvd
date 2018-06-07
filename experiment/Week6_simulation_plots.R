rm(list=ls())

# set up parameters
adj <- matrix(c(-.05, 0, -.05,
                0.15,-.1,-.1,
                -.25,-.15,0), 3, 3, byrow = T)
tmp <- svd(adj)
u_center <- t(tmp$u %*% diag(sqrt(tmp$d)))
v_center <- t(tmp$v %*% diag(sqrt(tmp$d)))

t(u_center) %*% v_center

u_num <- c(40, 20, 160)
u_label <- unlist(lapply(1:3, function(x){rep(x, u_num[x])}))
v_num <- c(60, 80, 200)
v_label <- unlist(lapply(1:3, function(x){rep(x, v_num[x])}))

# generate matrices
set.seed(10)
u_sig <- 0.05
u_dat <- do.call(rbind, lapply(1:3, function(x){
  MASS::mvrnorm(n = u_num[x], mu = u_center[,x], Sigma = u_sig*diag(3))
}))
v_sig <- 0.05
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
dropout_function <- dropout_create(.7,1)

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

length(which(setting == 1))/length(which(setting >= 1))

#############

.rotate <- function(mat){t(mat)[,nrow(mat):1]}

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

col_vec <- colorRamp_custom(c(0.584, 0.858, 0.564), c(0.803, 0.156, 0.211), 19)
col_vec <- c("white", col_vec)

.rotate <- function(mat){t(mat)[,nrow(mat):1]}

tmp <- as.numeric(dat)
tmp <- tmp[tmp!=0]

break_vec <- quantile(tmp, probs = seq(0, 1, length.out = 20))
break_vec <- c(-5, break_vec)

png("../figure/experiment/6_simulated_data.png", height = 1300, width = 2400, res = 300, units = "px")
par(mar = rep(0.5,4))
image(.rotate(dat), breaks = break_vec, col = col_vec, asp = nrow(dat)/ncol(dat),
      axes = F)

#put lines
row_idx <- 1-cumsum(u_num)[1:2]/n
col_idx <- cumsum(v_num)[1:2]/d

for(i in row_idx){
  lines(c(0,1), rep(i, 2), lwd = 6, col = "white")
  lines(c(0,1), rep(i, 2), lwd = 6, lty = 2)
}
for(i in col_idx){
  lines(rep(i, 2), c(0,1), lwd = 6, col = "white")
  lines(rep(i, 2), c(0,1), lwd = 6, lty = 2)
}

graphics.off()

################
res_svd <- svd(dat)
k <- 5
pred_dat <- res_svd$u[,1:k] %*% diag(res_svd$d[1:k]) %*% t(res_svd$v[,1:k])

zero_factor <- as.factor(as.numeric(dat != 0))
vec <- as.numeric(pred_dat)

naive_fit <- glm(zero_factor ~ vec, family=binomial(link='logit'))
naive_coef <- coef(naive_fit)

png("../figure/experiment/6_simulated_dropout_naive.png", height = 1300, width = 2400, res = 300, units = "px")
par(mfrow = c(1,2))

plot(as.numeric(mean_dat), as.numeric(pred_dat), xlab = "Mean value (uncensored)",
     ylab = "Predicted value", asp = T, col = rgb(0,0,0,0.2), pch = 16,
     main = "Latent mean value and\npredicted values")
lines(c(-100,100), rep(0,2), col = rgb(0.803, 0.156, 0.211), lwd = 2, lty = 2)
lines(rep(0,2), c(-100,100), col = rgb(0.803, 0.156, 0.211), lwd = 2, lty = 2)

tmp <- as.numeric(mean_dat)
spacing <- diff(range(tmp))/50
break_vec <- c(floor(min(tmp)/spacing):ceiling(max(tmp)/spacing))*spacing
hist(tmp[which(zero_factor == "0")], breaks = break_vec, col = rgb(0.803, 0.156, 0.211, 0.5),
     xlab = "Mean value (uncensored)", ylab = "Count", main = "Latent categorization of\nmean values")
hist(tmp[which(zero_factor == "1")], breaks = break_vec, col = rgb(0.584, 0.858, 0.564,0.5), add = T)
lines(rep(0,2), c(-100,1e10), lwd = 2, lty = 2)

graphics.off()

####################

set.seed(10)
fit_coef <- estimate_dropout(dat, threshold_quant_degree = 0.5, k = 3)
fit_coef

###

n <- nrow(dat); d <- ncol(dat)
adj_mat <- .form_adj_nonzero(dat)
g <- igraph::graph_from_adjacency_matrix(adj_mat, mode = "undirected")
igraph::V(g)$name <- as.character(1:(n+d))
idx <- as.numeric(.largest_average_degree_subgraph(g, threshold_quant = 0.5))
idx1 <- idx[which(idx <= n)]
idx2 <- idx[which(idx > n)] - n
pred_dat <- .funk_svd_prediction(dat[idx1, idx2], k = 3)
zero_factor <- as.factor(as.numeric(dat[idx1, idx2] != 0))
vec <- as.numeric(pred_dat)
tmp_mat <- data.frame(vec, zero_factor)
zz <- tmp_mat[(tmp_mat[,1] > 0), 1]
cutoff_val <- quantile(zz, probs = 0.5)

png("../figure/experiment/6_simulated_dropout_sophisticated.png", height = 1300, width = 2400, res = 300, units = "px")
par(mfrow = c(1,2))

plot(as.numeric(mean_dat[idx1, idx2]), as.numeric(pred_dat), xlab = "Mean value (uncensored)",
     ylab = "Predicted value", asp = T, col = rgb(0,0,0,0.2), pch = 16,
     main = "Latent mean value and\npredicted values")
lines(c(-100,100), rep(0,2), col = rgb(0.803, 0.156, 0.211), lwd = 2, lty = 2)
lines(rep(0,2), c(-100,100), col = rgb(0.803, 0.156, 0.211), lwd = 2, lty = 2)
lines(c(-100,100), rep(cutoff_val, 2), col = rgb(0.584, 0.858, 0.564), lwd = 2)

tmp <- as.numeric(mean_dat[idx1, idx2])
spacing <- diff(range(tmp))/50
break_vec <- c(floor(min(tmp)/spacing):ceiling(max(tmp)/spacing))*spacing
hist(tmp[which(zero_factor == "1")], breaks = break_vec, col = rgb(0.584, 0.858, 0.564,0.5),
     xlab = "Mean value (uncensored)", ylab = "Count", main = "Latent categorization of\nmean values")
hist(tmp[which(zero_factor == "0")], breaks = break_vec, col = rgb(0.803, 0.156, 0.211, 0.5), add = T)
lines(rep(0,2), c(-100,1e10), lwd = 2, lty = 2)

graphics.off()

######

tmp_svd <- svd(dat[idx1, idx2])
pred_dat_tmp <- tmp_svd$u[,1:3]%*%diag(tmp_svd$d[1:3])%*%t(tmp_svd$v[,1:3])
zero_factor <- as.factor(as.numeric(dat[idx1, idx2] != 0))
vec <- as.numeric(pred_dat_tmp)
tmp_mat <- data.frame(vec, zero_factor)
zz <- tmp_mat[(tmp_mat[,1] > 0), 1]
cutoff_val <- quantile(zz, probs = 0.5)
tmp_mat <- tmp_mat[which(tmp_mat[,1] > cutoff_val),]
fit <- stats::glm(zero_factor ~ vec, family=binomial(link='logit'), data = tmp_mat)
coef_vec <- stats::coef(fit)

######

x_seq <- seq(0, 4, length.out = 100)
true_prob <- sapply(x_seq, dropout_function)
naive_prob <- sapply(x_seq, dropout_create(naive_coef[1], naive_coef[2]))
good_prob <- sapply(x_seq, dropout_create(fit_coef[1], fit_coef[2]))

png("../figure/experiment/6_simulated_dropout_comparison.png", height = 1300, width = 1500, res = 300, units = "px")

plot(NA, xlab = "True value (unobserved)",
     ylab = "Probability of no dropout", main = "Comparison of dropout functions",
     xlim = c(0,4), ylim = c(0.5, 1))

lines(x_seq, true_prob, lwd = 5, lty = 2)
lines(x_seq, naive_prob, lwd = 5, col = rgb(0.803, 0.156, 0.211))
lines(x_seq, good_prob, lwd = 5, col = rgb(0.584, 0.858, 0.564))

legend("bottomright", c("True function", "Naive estimate", "Sophisticated estimate"),
       lty = 1, col = c("black",rgb(0.803, 0.156, 0.211), rgb(0.584, 0.858, 0.564)),
       lwd = 3)

graphics.off()

####################### # now for the weighted correlation matrices

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

weighted_gene_cor <- weighted_correlation(dat, fit_coef, 0.95)
gene_cor <- cor(dat)
true_cor <- cor(nodropout_dat)

col_idx <- sapply(1:(max(v_label)-1), function(x){
  length(which(v_label <= x))/length(v_label)
})

col_vec2 <- colorRampPalette(c(rgb(0.584, 0.858, 0.564), rgb(0.803, 0.156, 0.211)))(19)
break_vec <- quantile(c(as.numeric(true_cor), as.numeric(weighted_gene_cor),
                        as.numeric(gene_cor)), probs = seq(0, 1, length.out = 20))
.rotate <- function(mat){t(mat)[,nrow(mat):1]}

png(paste0("../figure/experiment/6_simulated_weighted_gene_correlation.png"), height = 2400, width = 7200, res = 300, units = "px")
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

sum((true_cor - gene_cor)^2)/nrow(dat)
sum((true_cor - weighted_gene_cor)^2)/nrow(dat)
sum((gene_cor - weighted_gene_cor)^2)/nrow(dat)

#### # try clustering

eig <- eigen(true_cor)
k <- 3
v_mat <- eig$vectors[,1:k] %*% diag(sqrt(eig$values[1:k]))
v_mat_spherical <- t(apply(v_mat, 1, function(x){x/.l2norm(x)}))

set.seed(10)
v_clust <- kmeans(v_mat_spherical, centers = 3, iter.max = 100, nstart = 10)
tab_true <- table(v_clust$cluster, v_label)
##
eig <- eigen(weighted_gene_cor)
k <- 3
v_mat <- eig$vectors[,1:k] %*% diag(sqrt(eig$values[1:k]))
v_mat_spherical <- t(apply(v_mat, 1, function(x){x/.l2norm(x)}))

set.seed(10)
v_clust <- kmeans(v_mat_spherical, centers = 3, iter.max = 100, nstart = 10)
tab_weighted <- table(v_clust$cluster, v_label)

########### # do the same but for cells

weighted_cell_cor <- weighted_correlation(t(dat), fit_coef, 0.95)
cell_cor <- cor(t(dat))
true_cell_cor <- cor(t(nodropout_dat))

col_idx <- sapply(1:(max(u_label)-1), function(x){
  length(which(u_label <= x))/length(u_label)
})

col_vec2 <- colorRampPalette(c(rgb(0.584, 0.858, 0.564), rgb(0.803, 0.156, 0.211)))(19)
break_vec <- quantile(c(as.numeric(true_cell_cor), as.numeric(weighted_cell_cor),
                        as.numeric(cell_cor)), probs = seq(0, 1, length.out = 20))
.rotate <- function(mat){t(mat)[,nrow(mat):1]}

png(paste0("../figure/experiment/6_simulated_weighted_cell_correlation.png"), height = 2400, width = 7200, res = 300, units = "px")
par(mfrow = c(1,3))
image(.rotate(true_cell_cor), breaks = break_vec, col = col_vec2, asp = T, axes = F,
      main = "True correlation", cex.main = 3)

for(j in col_idx){
  lines(c(0,1), rep(1-j, 2), lwd = 6, col = "white")
  lines(c(0,1), rep(1-j, 2), lwd = 6, lty = 2)
}
for(j in col_idx){
  lines(rep(j, 2), c(0,1), lwd = 6, col = "white")
  lines(rep(j, 2), c(0,1), lwd = 6, lty = 2)
}

image(.rotate(cell_cor), breaks = break_vec, col = col_vec2, asp = T, axes = F,
      main = "Estimated correlation", cex.main = 3)

for(j in col_idx){
  lines(c(0,1), rep(1-j, 2), lwd = 6, col = "white")
  lines(c(0,1), rep(1-j, 2), lwd = 6, lty = 2)
}
for(j in col_idx){
  lines(rep(j, 2), c(0,1), lwd = 6, col = "white")
  lines(rep(j, 2), c(0,1), lwd = 6, lty = 2)
}

image(.rotate(weighted_cell_cor), breaks = break_vec, col = col_vec2, asp = T, axes = F,
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

sum((true_cor - gene_cor)^2)/nrow(dat)
sum((true_cor - weighted_gene_cor)^2)/nrow(dat)
sum((gene_cor - weighted_gene_cor)^2)/nrow(dat)

#### # try clustering

eig <- eigen(true_cell_cor)
k <- 3
u_mat <- eig$vectors[,1:k] %*% diag(sqrt(eig$values[1:k]))
u_mat_spherical <- t(apply(u_mat, 1, function(x){x/.l2norm(x)}))

set.seed(10)
u_clust <- kmeans(u_mat_spherical, centers = 3, iter.max = 100, nstart = 10)
tab_true <- table(u_clust$cluster, u_label)
##
eig <- eigen(weighted_cell_cor)
k <- 3
u_mat <- eig$vectors[,1:k] %*% diag(sqrt(eig$values[1:k]))
u_mat_spherical <- t(apply(u_mat, 1, function(x){x/.l2norm(x)}))

set.seed(10)
u_clust <- kmeans(u_mat_spherical, centers = 3, iter.max = 100, nstart = 10)
tab_weighted <- table(u_clust$cluster, u_label)
