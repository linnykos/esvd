rm(list=ls())
load("../../SOUP/data/camp.rda")

dat <- camp$counts

idx <- which(colnames(dat) %in% camp$select.genes)
dat <- dat[,idx]

dat <- t(apply(dat, 1, function(x){x/sum(x)}))
dat <- 10*log(dat + 1)

dim(dat)

##################

res_svd <- svd(dat)

plot(res_svd$d[1:50])
k <- 4
u_mat <- res_svd$u[,1:k] %*% diag(sqrt(res_svd$d[1:k]))
v_mat <- res_svd$v[,1:k] %*% diag(sqrt(res_svd$d[1:k]))

plot(u_mat[,1], u_mat[,2], asp = T, xlim = range(c(u_mat[,1], 0)), ylim = range(c(u_mat[,2], 0)))
plot(v_mat[,1], v_mat[,2], asp = T, xlim = range(c(v_mat[,1], 0)), ylim = range(c(v_mat[,2], 0)))

u_mat_spherical <- t(apply(u_mat, 1, function(x){x/.l2norm(x)}))
v_mat_spherical <- t(apply(v_mat, 1, function(x){x/.l2norm(x)}))

set.seed(10)
u_clust <- kmeans(u_mat_spherical, centers = 3, iter.max = 100, nstart = 10)
v_clust <- kmeans(v_mat_spherical, centers = 3, iter.max = 100, nstart = 10)

#####################


cell_type <- as.numeric(as.factor(camp$cell.info[,2]))
row_idx <- unlist(lapply(1:max(cell_type), function(x){
  which(cell_type == x)
}))
col_idx <- unlist(lapply(1:max(u_clust$cluster), function(x){
  sort(which(v_clust$cluster == x))
}))

dat3 <- dat[row_idx, col_idx]

tmp <- as.numeric(dat3)
tmp <- tmp[tmp!=0]

break_vec <- quantile(tmp, probs = seq(0, 1, length.out = 20))
break_vec <- c(-5, break_vec)

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

cell_type <- camp$cell.info[,2]
cell_type_coarse <- as.numeric(as.factor(sapply(cell_type, function(x){substr(x, 1, 1)})))
row_idx <- unlist(lapply(1:max(cell_type_coarse), function(x){
  which(cell_type_coarse == x)
}))

dat4 <- dat[row_idx, col_idx]

png("../figure/experiment/6_camp_data_sorted_3group.png", height = 1300, width = 2400, res = 300, units = "px")
par(mar = rep(0.5,4))
image(.rotate(dat4), breaks = break_vec, col = col_vec, asp = nrow(dat4)/ncol(dat4),
      axes = F)

#put lines
row_idx_lines <- sapply(1:(max(cell_type_coarse)-1), function(x){
  1-length(which(cell_type_coarse <= x))/length(cell_type_coarse)
})
col_idx_lines <- sapply(1:(max(v_clust$cluster)-1), function(x){
  length(which(v_clust$cluster <= x))/length(v_clust$cluster)
})

for(i in row_idx_lines){
  lines(c(0,1), rep(i, 2), lwd = 6, col = "white")
  lines(c(0,1), rep(i, 2), lwd = 6, lty = 2)
}
for(i in col_idx_lines){
  lines(rep(i, 2), c(0,1), lwd = 6, col = "white")
  lines(rep(i, 2), c(0,1), lwd = 6, lty = 2)
}

graphics.off()

############

set.seed(10)
fit_coef <- estimate_dropout(dat, threshold_quant_degree = 0.5, k = 4,
                             threshold_quant_logistic = 0.1)
fit_coef

x_seq <- seq(0, 1.5, length.out = 100)
good_prob <- sapply(x_seq, dropout_create(fit_coef[1], fit_coef[2]))

png("../figure/experiment/6_camp_dropout.png", height = 1300, width = 2400, res = 300, units = "px")
par(mfrow = c(1,2))
plot(NA, xlab = "True value (unobserved)",
     ylab = "Probability of no dropout", main = "Comparison of dropout functions",
     xlim = c(0,1.5), ylim = c(0.4, 1))

lines(x_seq, good_prob, lwd = 5, col = rgb(0.584, 0.858, 0.564))

hist(as.numeric(dat), breaks = 50, col = "gray", ylim = c(0,4000), xlim = c(0, 1.5),
     main = "Histogram of values", xlab = "Values")

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

weighted_cell_cor <- weighted_correlation(t(dat4), fit_coef, 0.95)
cell_cor <- cor(t(dat4))

col_idx <- sapply(1:(max(cell_type_coarse)-1), function(x){
  length(which(cell_type_coarse <= x))/length(cell_type_coarse)
})

col_vec2 <- colorRampPalette(c(rgb(0.584, 0.858, 0.564), rgb(0.803, 0.156, 0.211)))(19)
break_vec <- quantile(c(as.numeric(weighted_cell_cor),
                        as.numeric(cell_cor)), probs = seq(0, 1, length.out = 20))
.rotate <- function(mat){t(mat)[,nrow(mat):1]}

png(paste0("../figure/experiment/6_camp_weighted_cell_correlation.png"), height = 2400, width = 4800, res = 300, units = "px")
par(mfrow = c(1,2))

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

###################

# try clustering

eig <- eigen(cell_cor)
k <- 3
u_mat <- eig$vectors[,1:k] %*% diag(sqrt(eig$values[1:k]))
u_mat_spherical <- t(apply(u_mat, 1, function(x){x/.l2norm(x)}))

set.seed(10)
u_clust <- kmeans(u_mat_spherical, centers = 3, iter.max = 100, nstart = 10)
tab_true <- table(cell_type_coarse, u_clust$cluster)
##
eig <- eigen(weighted_cell_cor)
k <- 3
u_mat <- eig$vectors[,1:k] %*% diag(sqrt(eig$values[1:k]))
u_mat_spherical <- t(apply(u_mat, 1, function(x){x/.l2norm(x)}))

set.seed(10)
u_clust <- kmeans(u_mat_spherical, centers = 3, iter.max = 100, nstart = 10)
tab_weighted <- table(cell_type_coarse, u_clust$cluster)
