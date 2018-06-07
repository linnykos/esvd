rm(list=ls())

# set up parameters
adj <- matrix(c(-.1, 0, -.1,
                0.15,-.25,-.25,
                -.75,-.25,0), 3, 3, byrow = T)
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

dropout_create <- function(a, b){
  function(x){ 1/(1+exp(-(a+b*x))) }
}
dropout_function <- dropout_create(.7,1)

plot(seq(-5,5,length.out = 100), dropout_function(seq(-5,5,length.out = 100)))
lines(rep(0,2), c(-100,100), col = "red", lwd = 2)

set.seed(10)
for(i in 1:n){
  for(j in 1:d){
    if(mean_dat[i,j] <= 0) {
      dat[i,j] <- 0
      setting[i,j] <- 0
    } else {
      val <- rexp(1, rate = 1/mean_dat[i,j])
      bool <- rbinom(1, 1, prob = dropout_function(val))
      setting[i,j] <- bool + 1
      dat[i,j] <- bool*val
    }
  }
}

length(which(setting == 1))/length(which(setting >= 1))


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

############################

# let's try using a standard SVD to make a crude estimate of the dropout function
res_svd <- svd(dat)
k <- 5
pred_dat <- res_svd$u[,1:k] %*% diag(res_svd$d[1:k]) %*% t(res_svd$v[,1:k])
plot(as.numeric(mean_dat), as.numeric(pred_dat), col = rgb(0,0,0,0.2))

zero_factor <- as.factor(as.numeric(dat != 0))
vec <- as.numeric(pred_dat)

fit <- glm(zero_factor ~ vec, family=binomial(link='logit'))

actual_prob <- dropout_function(as.numeric(mean_dat))
fit_prob <- dropout_create(coef(fit)[1], coef(fit)[2])(as.numeric(mean_dat))

plot(actual_prob, fit_prob, col = rgb(0,0,0,0.2), asp = T, xlim = c(0,1), ylim = c(0,1))
lines(c(0,1), c(0,1), col = "red", lwd = 2)
sum((actual_prob - fit_prob)^2)/length(actual_prob)

### what if we tried to fix the "zero_factor" by clustering and seeing which blocks are 0?
u_mat <- res_svd$u[,1:k] %*% diag(sqrt(res_svd$d[1:k]))
v_mat <- res_svd$v[,1:k] %*% diag(sqrt(res_svd$d[1:k]))

u_mat_spherical <- t(apply(u_mat, 1, function(x){x/.l2norm(x)}))
v_mat_spherical <- t(apply(v_mat, 1, function(x){x/.l2norm(x)}))

set.seed(10)
u_clust <- kmeans(u_mat_spherical, centers = 3, iter.max = 100, nstart = 10)
v_clust <- kmeans(v_mat_spherical, centers = 3, iter.max = 100, nstart = 10)

table(u_label, u_clust$cluster)
table(v_label, v_clust$cluster)

for(i in 1:3){
  for(j in 1:3){
    idx1 <- which(u_clust$cluster == i)
    idx2 <- which(v_clust$cluster == j)
    print(paste0(i, ", ", j, ": ", length(which(dat[idx1, idx2] != 0))/(length(idx1)*length(idx2))))
  }
}

pred_setting <- matrix(0, n, d)

for(i in 1:3){
  for(j in 1:3){
    idx1 <- which(u_clust$cluster == i)
    idx2 <- which(v_clust$cluster == j)
    percentage <- length(which(dat[idx1, idx2] != 0))/(length(idx1)*length(idx2))

    if(percentage >= 0.2) pred_setting[idx1, idx2] <- 1 else pred_setting[idx1, idx2] <- 0
  }
}

zero_factor <- as.factor(as.numeric(pred_setting != 0))
vec <- as.numeric(pred_dat)

fit <- glm(zero_factor ~ vec, family=binomial(link='logit'))

actual_prob <- dropout_function(as.numeric(mean_dat))
fit_prob <- dropout_create(coef(fit)[1], coef(fit)[2])(as.numeric(mean_dat))

plot(actual_prob, fit_prob, col = rgb(0,0,0,0.2), asp = T, xlim = c(0,1), ylim = c(0,1))
lines(c(0,1), c(0,1), col = "red", lwd = 2)
sum((actual_prob - fit_prob)^2)/length(actual_prob)


#############################

# attempt number 2: matrix factorization
dat2 <- dat
dat2[dat2 == 0] <- NA
dat2 <- as(dat2, "realRatingMatrix")
dat2
set.seed(10)
funk_svd <- recommenderlab::funkSVD(dat2, k = 5)
pred_dat <- tcrossprod(funk_svd$U, funk_svd$V)

plot(as.numeric(mean_dat), as.numeric(pred_dat), col = rgb(0,0,0,0.2))

zero_factor <- as.factor(as.numeric(dat != 0))
vec <- as.numeric(pred_dat)

fit <- glm(zero_factor ~ vec, family=binomial(link='logit'))

actual_prob <- dropout_function(as.numeric(mean_dat))
fit_prob <- dropout_create(coef(fit)[1], coef(fit)[2])(as.numeric(mean_dat))

plot(actual_prob, fit_prob, col = rgb(0,0,0,0.2))
sum((actual_prob - fit_prob)^2)/length(actual_prob)

# play with different settings, to get some negatives
funk_svd <- recommenderlab::funkSVD(dat2, k = 5, lambda = 0.01)
pred_dat <- tcrossprod(funk_svd$U, funk_svd$V)

plot(as.numeric(mean_dat), as.numeric(pred_dat), col = rgb(0,0,0,0.2))

zero_factor <- as.factor(as.numeric(dat != 0))
vec <- as.numeric(pred_dat)

fit <- glm(zero_factor ~ vec, family=binomial(link='logit'))

actual_prob <- dropout_function(as.numeric(mean_dat))
fit_prob <- dropout_create(coef(fit)[1], coef(fit)[2])(as.numeric(mean_dat))

plot(actual_prob, fit_prob, col = rgb(0,0,0,0.2), asp = T, xlim = c(0,1), ylim = c(0,1))
lines(c(0,1), c(0,1), col = "red", lwd = 2)
sum((actual_prob - fit_prob)^2)/length(actual_prob)

#####################

# let's see how good it would be if mean_dat were given
dat_neg <- matrix(0, n, d)

set.seed(10)
for(i in 1:n){
  for(j in 1:d){
    bool <- rbinom(1, 1, prob = dropout_function(mean_dat[i,j]))
    dat_neg[i,j] <- bool*mean_dat[i,j]
  }
}

res_svd <- svd(dat_neg)
k <- 5
pred_dat <- res_svd$u[,1:k] %*% diag(res_svd$d[1:k]) %*% t(res_svd$v[,1:k])
plot(as.numeric(mean_dat), as.numeric(pred_dat), col = rgb(0,0,0,0.2), asp = T)

zero_factor <- as.factor(as.numeric(dat_neg != 0))
vec <- as.numeric(pred_dat)

fit <- glm(zero_factor ~ vec, family=binomial(link='logit'))

actual_prob <- dropout_function(as.numeric(mean_dat))
fit_prob <- dropout_create(coef(fit)[1], coef(fit)[2])(as.numeric(mean_dat))

plot(actual_prob, fit_prob, col = rgb(0,0,0,0.2), asp = T, xlim = c(0,1), ylim = c(0,1))
lines(c(0,1), c(0,1), col = "red", lwd = 2)
sum((actual_prob - fit_prob)^2)/length(actual_prob)

#############

# let's try to see if we can estimate the probability from just one block?

for(i in 1:3){
  for(j in 1:3){
    idx1 <- which(u_clust$cluster == i)
    idx2 <- which(v_clust$cluster == j)
    percentage <- length(which(dat[idx1, idx2] != 0))/(length(idx1)*length(idx2))

    if(percentage >= 0.2) {
      dat2 <- dat[idx1, idx2]
      dat2[dat2 == 0] <- NA
      dat2 <- as(dat2, "realRatingMatrix")

      funk_svd <- recommenderlab::funkSVD(dat2, k = 5, lambda = 0.01)
      pred_dat <- tcrossprod(funk_svd$U, funk_svd$V)

      zero_factor <- as.factor(as.numeric(dat[idx1, idx2] != 0))
      vec <- as.numeric(pred_dat)

      fit <- glm(zero_factor ~ vec, family=binomial(link='logit'))

      actual_prob <- dropout_function(as.numeric(mean_dat[idx1, idx2]))
      fit_prob <- dropout_create(coef(fit)[1], coef(fit)[2])(as.numeric(mean_dat[idx1, idx2]))

      plot(actual_prob, fit_prob, col = rgb(0,0,0,0.2), asp = T, xlim = c(0,1), ylim = c(0,1))
      lines(c(0,1), c(0,1), col = "red", lwd = 2)
      print(paste0(i, ", ", j, ": ", round(sum((actual_prob - fit_prob)^2)/length(actual_prob),5)))
    }
  }
}

# zoom in for some analysis
i <- 3; j <- 2
idx1 <- which(u_clust$cluster == i)
idx2 <- which(v_clust$cluster == j)
table(setting[idx1, idx2])
plot(mean_dat[idx1, idx2], zero_factor)

tmp <- as.numeric(mean_dat[idx1, idx2])
break_vec <- seq(min(tmp), max(tmp), length.out = 50)
hist(tmp[which(zero_factor == "0")], breaks = break_vec, col = rgb(1,0,0,0.2))
hist(tmp[which(zero_factor == "1")], breaks = break_vec, col = rgb(0,0,1,0.2), add = T)

######### #does knowing exactly the correct block help at all? ... doesn't seem like it

idx1 <- which(u_label == 3)
idx2 <- which(v_label == 3)
dat2 <- dat[idx1, idx2]
dat2[dat2 == 0] <- NA
dat2 <- as(dat2, "realRatingMatrix")

funk_svd <- recommenderlab::funkSVD(dat2, k = 5, lambda = 0.01)
pred_dat <- tcrossprod(funk_svd$U, funk_svd$V)

zero_factor <- as.factor(as.numeric(dat[idx1, idx2] != 0))
vec <- as.numeric(pred_dat)

fit <- glm(zero_factor ~ vec, family=binomial(link='logit'))

actual_prob <- dropout_function(as.numeric(mean_dat[idx1, idx2]))
fit_prob <- dropout_create(coef(fit)[1], coef(fit)[2])(as.numeric(mean_dat[idx1, idx2]))

plot(actual_prob, fit_prob, col = rgb(0,0,0,0.2), asp = T, xlim = c(0,1), ylim = c(0,1))
lines(c(0,1), c(0,1), col = "red", lwd = 2)

table(zero_factor)
table(as.numeric(setting[idx1, idx2]))

######## # let's instead find an oracle that has the /least/ zero in the setting
# let's find the densest component of the bipartiate graph where setting of 1 or 2 is an edge
adj_mat <- matrix(0, n+d, n+d)
for(i in 1:n){
  for(j in 1:d){
    adj_mat[i, n+j] <- ifelse(setting[i,j] == 0, 0, 1)
    adj_mat[n+j, i] <- adj_mat[i, n+j]
  }
}

library(igraph)
g <- igraph::graph_from_adjacency_matrix(adj_mat, mode = "undirected")
igraph::V(g)$name <- as.character(1:(n+d))

# now find the subgraph with the highest degree
vertex_list <- vector("list", n+d)
num_edges_list <- rep(0, n+d)

set.seed(10)
g_tmp <- g
for(i in 1:length(vertex_list)){
  name_vec <- igraph::V(g_tmp)$name
  vec <- igraph::degree(g_tmp)
  min_val <- min(vec)
  idx <- which(vec == min_val)
  if(length(idx) > 1) idx <- idx[sample(length(idx), 1)]

  name_vec <- name_vec[-idx]
  g_tmp <- igraph::induced_subgraph(g_tmp, name_vec)

  vertex_list[[i]] <- name_vec
  num_edges_list[i] <- igraph::ecount(g_tmp)/igraph::vcount(g_tmp)
}

plot(num_edges_list)

idx <- which.max(num_edges_list)
vertex_vec <- as.numeric(vertex_list[[idx]])
idx1 <- vertex_vec[which(vertex_vec <= n)]
idx2 <- vertex_vec[which(vertex_vec > n)] - n

table(as.numeric(setting[idx1, idx2]))

## funk svd

dat2 <- dat[idx1, idx2]
dat2[dat2 == 0] <- NA
dat2 <- as(dat2, "realRatingMatrix")

funk_svd <- recommenderlab::funkSVD(dat2, k = 5, lambda = 0.01)
pred_dat <- tcrossprod(funk_svd$U, funk_svd$V)

plot(as.numeric(mean_dat[idx1, idx2]), as.numeric(pred_dat))

zero_factor <- as.factor(as.numeric(dat[idx1, idx2] != 0))
vec <- as.numeric(pred_dat)

fit <- glm(zero_factor ~ vec, family=binomial(link='logit'))

actual_prob <- dropout_function(as.numeric(mean_dat[idx1, idx2]))
fit_prob <- dropout_create(coef(fit)[1], coef(fit)[2])(as.numeric(mean_dat[idx1, idx2]))

plot(actual_prob, fit_prob, col = rgb(0,0,0,0.2), asp = T, xlim = c(0,1), ylim = c(0,1))
lines(c(0,1), c(0,1), col = "red", lwd = 2)

tmp <- as.numeric(mean_dat[idx1, idx2])
break_vec <- seq(min(tmp), max(tmp), length.out = 50)
hist(tmp[which(zero_factor == "0")], breaks = break_vec, col = rgb(1,0,0,0.2))
hist(tmp[which(zero_factor == "1")], breaks = break_vec, col = rgb(0,0,1,0.2), add = T)

# what about just using svd?
res_svd_tmp <-svd(dat[idx1, idx2])
pred_dat <- res_svd_tmp$u[,1:k] %*% diag(res_svd_tmp$d[1:k]) %*% t(res_svd_tmp$v[,1:k])

zero_factor <- as.factor(as.numeric(dat[idx1, idx2] != 0))
vec <- as.numeric(pred_dat)

fit <- glm(zero_factor ~ vec, family=binomial(link='logit'))

actual_prob <- dropout_function(as.numeric(mean_dat[idx1, idx2]))
fit_prob <- dropout_create(coef(fit)[1], coef(fit)[2])(as.numeric(mean_dat[idx1, idx2]))

plot(actual_prob, fit_prob, col = rgb(0,0,0,0.2), asp = T, xlim = c(0,1), ylim = c(0,1))
lines(c(0,1), c(0,1), col = "red", lwd = 2)

# try dropping the negative ones for predicting the dropout function
plot(as.numeric(mean_dat[idx1, idx2]), as.numeric(pred_dat))
tmp_mat <- data.frame(vec, zero_factor)
tmp_mat <- tmp_mat[which(tmp_mat[,1] > 0.15),]
fit <- glm(zero_factor ~ vec, family=binomial(link='logit'), data = tmp_mat)
coef(fit)

actual_prob <- dropout_function(seq(-5,5,length.out = 100))
fit_prob <- dropout_create(coef(fit)[1], coef(fit)[2])(seq(-5,5,length.out = 100))
idx <- which.min(abs(seq(-5,5,length.out = 100)))

col_vec <- rep(rgb(0,0,0,0.2), length(fit_prob))
col_vec[idx] <- "red"
plot(actual_prob, fit_prob, col = col_vec, asp = T, xlim = c(0,1), ylim = c(0,1))
lines(c(0,1), c(0,1), col = "red", lwd = 2)

plot(tmp_mat[,1], predict.glm(fit, newdata = tmp_mat, type = "response"))
plot(tmp_mat[,1], sapply(tmp_mat[,1], dropout_function))

###### # combining different ideas: funk_svd + threshold. This looks great

dat2 <- dat[idx1, idx2]
dat2[dat2 == 0] <- NA
dat2 <- as(dat2, "realRatingMatrix")

funk_svd <- recommenderlab::funkSVD(dat2, k = 5, lambda = 0.01)
pred_dat <- tcrossprod(funk_svd$U, funk_svd$V)

zero_factor <- as.factor(as.numeric(dat[idx1, idx2] != 0))
vec <- as.numeric(pred_dat)

tmp_mat <- data.frame(vec, zero_factor)
zz <- tmp_mat[(tmp_mat[,1] > 0),1]
cutoff_val <- quantile(zz, probs = 0.5)
tmp_mat <- tmp_mat[which(tmp_mat[,1] > cutoff_val),]
fit <- glm(zero_factor ~ vec, family=binomial(link='logit'), data = tmp_mat)
coef(fit)

actual_prob <- dropout_function(seq(-5,5,length.out = 100))
fit_prob <- dropout_create(coef(fit)[1], coef(fit)[2])(seq(-5,5,length.out = 100))
idx <- which.min(abs(seq(-5,5,length.out = 100)))

col_vec <- rep(rgb(0,0,0,0.2), length(fit_prob))
col_vec[idx] <- "red"
plot(actual_prob, fit_prob, col = col_vec, asp = T, xlim = c(0,1), ylim = c(0,1))
lines(c(0,1), c(0,1), col = "red", lwd = 2)

plot(tmp_mat[,1], predict.glm(fit, newdata = tmp_mat, type = "response"))
plot(tmp_mat[,1], sapply(tmp_mat[,1], dropout_function))

###################

# now try it on the data:
estimate_dropout_parameters <- function(dat){
  n <- nrow(dat); d <- ncol(dat)

  adj_mat <- matrix(0, n+d, n+d)
  for(i in 1:n){
    for(j in 1:d){
      adj_mat[i, n+j] <- ifelse(dat[i,j] == 0, 0, 1)
      adj_mat[n+j, i] <- adj_mat[i, n+j]
    }
  }

  library(igraph)
  g <- igraph::graph_from_adjacency_matrix(adj_mat, mode = "undirected")
  igraph::V(g)$name <- as.character(1:(n+d))

  # now find the subgraph with the highest degree
  vertex_list <- vector("list", n+d)
  num_edges_list <- rep(0, n+d)

  set.seed(10)
  g_tmp <- g
  for(i in 1:length(vertex_list)){
    name_vec <- igraph::V(g_tmp)$name
    vec <- igraph::degree(g_tmp)
    min_val <- min(vec)
    idx <- which(vec == min_val)
    if(length(idx) > 1) idx <- idx[sample(length(idx), 1)]

    name_vec <- name_vec[-idx]
    g_tmp <- igraph::induced_subgraph(g_tmp, name_vec)

    vertex_list[[i]] <- name_vec
    num_edges_list[i] <- igraph::ecount(g_tmp)/igraph::vcount(g_tmp)
  }

  plot(num_edges_list)

  idx <- which.max(num_edges_list)
  vertex_vec <- as.numeric(vertex_list[[idx]])
  idx1 <- vertex_vec[which(vertex_vec <= n)]
  idx2 <- vertex_vec[which(vertex_vec > n)] - n

  dat2 <- dat[idx1, idx2]
  dat2[dat2 == 0] <- NA
  dat2 <- as(dat2, "realRatingMatrix")

  funk_svd <- recommenderlab::funkSVD(dat2, k = 5, lambda = 0.01)
  pred_dat <- tcrossprod(funk_svd$U, funk_svd$V)

  zero_factor <- as.factor(as.numeric(dat[idx1, idx2] != 0))
  vec <- as.numeric(pred_dat)

  tmp_mat <- data.frame(vec, zero_factor)
  zz <- tmp_mat[(tmp_mat[,1] > 0),1]
  cutoff_val <- quantile(zz, probs = 0.5)
  tmp_mat <- tmp_mat[which(tmp_mat[,1] > cutoff_val),]
  fit <- glm(zero_factor ~ vec, family=binomial(link='logit'), data = tmp_mat)
  coef(fit)
}

fit_coef <- estimate_dropout_parameters(dat)

actual_prob <- dropout_function(seq(-5,5,length.out = 100))
fit_prob <- dropout_create(fit_coef[1], fit_coef[2])(seq(-5,5,length.out = 100))
idx <- which.min(abs(seq(-5,5,length.out = 100)))

col_vec <- rep(rgb(0,0,0,0.2), length(fit_prob))
col_vec[idx] <- "red"
plot(actual_prob, fit_prob, col = col_vec, asp = T, xlim = c(0,1), ylim = c(0,1))
lines(c(0,1), c(0,1), col = "red", lwd = 2)

actual_prob <- dropout_function(seq(0,5,length.out = 100))
fit_prob <- dropout_create(fit_coef[1], fit_coef[2])(seq(0,5,length.out = 100))
plot(seq(0,5,length.out=100), actual_prob, pch = 16, ylim = range(c(actual_prob, fit_prob)))
points(seq(0,5,length.out=100), fit_prob, col = "red", pch = 16)

###########################

# now let's see about this weighted covariance matrix
library(boot)

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

# also try a different (the wrong) dropout coefficients
res_svd <- svd(dat)
k <- 5
pred_dat <- res_svd$u[,1:k] %*% diag(res_svd$d[1:k]) %*% t(res_svd$v[,1:k])

zero_factor <- as.factor(as.numeric(dat != 0))
vec <- as.numeric(pred_dat)

fit2 <- glm(zero_factor ~ vec, family=binomial(link='logit'))
fit2_coef <- coef(fit2)
weighted_gene_cor2 <- weighted_correlation(dat, fit2_coef, 0.95)

# plot all three

col_idx <- sapply(1:(max(v_label)-1), function(x){
  length(which(v_label <= x))/length(v_label)
})

col_vec2 <- colorRampPalette(c(rgb(0.584, 0.858, 0.564), rgb(0.803, 0.156, 0.211)))(19)
break_vec <- quantile(c(as.numeric(weighted_gene_cor2), as.numeric(weighted_gene_cor),
                        as.numeric(gene_cor)), probs = seq(0, 1, length.out = 20))
.rotate <- function(mat){t(mat)[,nrow(mat):1]}

png(paste0("../figure/experiment/6_simulated_weighted_gene_correlation.png"), height = 2400, width = 7200, res = 300, units = "px")
par(mfrow = c(1,3))
image(.rotate(gene_cor), breaks = break_vec, col = col_vec2, asp = T, axes = F)

for(j in col_idx){
  lines(c(0,1), rep(1-j, 2), lwd = 6, col = "white")
  lines(c(0,1), rep(1-j, 2), lwd = 6, lty = 2)
}
for(j in col_idx){
  lines(rep(j, 2), c(0,1), lwd = 6, col = "white")
  lines(rep(j, 2), c(0,1), lwd = 6, lty = 2)
}

image(.rotate(weighted_gene_cor), breaks = break_vec, col = col_vec2, asp = T, axes = F)

for(j in col_idx){
  lines(c(0,1), rep(1-j, 2), lwd = 6, col = "white")
  lines(c(0,1), rep(1-j, 2), lwd = 6, lty = 2)
}
for(j in col_idx){
  lines(rep(j, 2), c(0,1), lwd = 6, col = "white")
  lines(rep(j, 2), c(0,1), lwd = 6, lty = 2)
}

image(.rotate(weighted_gene_cor2), breaks = break_vec, col = col_vec2, asp = T, axes = F)

for(j in col_idx){
  lines(c(0,1), rep(1-j, 2), lwd = 6, col = "white")
  lines(c(0,1), rep(1-j, 2), lwd = 6, lty = 2)
}
for(j in col_idx){
  lines(rep(j, 2), c(0,1), lwd = 6, col = "white")
  lines(rep(j, 2), c(0,1), lwd = 6, lty = 2)
}

graphics.off()

sum((gene_cor - weighted_gene_cor)^2)/nrow(dat)
sum((gene_cor - weighted_gene_cor2)^2)/nrow(dat)

eig <- eigen(gene_cor)
k <- 5
v_mat <- eig$vectors[,1:k] %*% diag(sqrt(eig$values[1:k]))
v_mat_spherical <- t(apply(v_mat, 1, function(x){x/.l2norm(x)}))

set.seed(10)
v_clust <- kmeans(v_mat_spherical, centers = 3, iter.max = 100, nstart = 10)
table(v_clust$cluster, v_label)

###########################


