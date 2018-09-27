rm(list=ls())
load("../../SOUP/data/camp.rda")

dat <- camp$counts

idx <- which(colnames(dat) %in% camp$select.genes)
dat <- dat[,idx]

dat <- t(apply(dat, 1, function(x){x/sum(x)}))
dat <- 10*log(dat + 1)

dim(dat)
dat <- dat[order(camp$cell.info[,2]),]
camp$cell.info <- camp$cell.info[order(camp$cell.info[,2]),]

length(which(dat == 0))/prod(dim(dat))

#####################

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

cell_type <- camp$cell.info[,2]
cell_type_coarse <- as.numeric(as.factor(sapply(cell_type, function(x){substr(x, 1, 1)})))
row_idx <- unlist(lapply(1:max(cell_type_coarse), function(x){
  which(cell_type_coarse == x)
}))
row_idx_lines <- sapply(1:(max(cell_type_coarse)-1), function(x){
  1-length(which(cell_type_coarse <= x))/length(cell_type_coarse)
})

png("../figure/experiment/21_camp_data.png", height = 800, width = 1500, res = 300, units = "px")
par(mar = rep(0.5,4))
image(.rotate(dat), breaks = break_vec, col = col_vec, asp = nrow(dat)/ncol(dat),
      axes = F)

for(i in row_idx_lines){
  lines(c(0,1), rep(i, 2), lwd = 3, col = "white")
  lines(c(0,1), rep(i, 2), lwd = 3, lty = 2)
}
graphics.off()

##################

dropout_res <- vector("list", length = ncol(dat))
for(i in 1:ncol(dat)){
  set.seed(10)
  if(i %% floor(ncol(dat)/10) == 0) cat('*')
  x <- dat[,i]
  if(all(x == 0)) {
    dropout_res[[i]] <- NA
    next()
  }
  if(length(x[x > 0]) <= 10) {
    dropout_res[[i]] <- NA
    next()
  }

  dropout_res[[i]] <- tryCatch({
    .em_mixture(x)
  }, error = function(e){
    print(paste("Error in ", i))
    return(NA)
  })
}

## a function to weigh each element of the matrix
## weight of 0 means it is VERY likely a dropout
weight_mat <- matrix(NA, ncol = ncol(dat), nrow = nrow(dat))
for(i in 1:ncol(weight_mat)){
  if(any(is.na(dropout_res[[i]]))){
    weight_mat[,i] <- 0
  } else {
    set.seed(10)
    x <- dat[,i]
    x2 <- .jitter_zeros(x)
    wt <- .calculate_weight(x2, dropout_res[[i]])
    weight_mat[,i] <- wt[,2]/(wt[,1] + wt[,2])
  }
}

col_vec <- colorRampPalette(c("blue3", "gold"))(20)
png("../figure/experiment/21_camp_weight.png", height = 800, width = 1500, res = 300, units = "px")
par(mar = rep(0.5,4))
image(.rotate(weight_mat), col = col_vec, asp = nrow(weight_mat)/ncol(weight_mat),
      axes = F)

for(i in row_idx_lines){
  lines(c(0,1), rep(i, 2), lwd = 3, col = "white")
  lines(c(0,1), rep(i, 2), lwd = 3, lty = 2)
}
graphics.off()

dropout_idx <- vector("list", length = ncol(dat))
for(i in 1:ncol(dat)){
  if(!any(is.na(dropout_res[[i]]))){
    set.seed(10)
    x <- dat[,i]
    x2 <- .jitter_zeros(x)
    dropout_idx[[i]] <- which(.compute_dropout(dropout_res[[i]], x2) == 1)
  }
}

#1 means keep, 0 means drop
dropout_mat <- matrix(1, ncol = ncol(dat), nrow = nrow(dat))
for(i in 1:ncol(dat)){
  if(any(is.na(dropout_res[[i]]))){
    dropout_mat[,i] <- 0
  } else {
    dropout_mat[dropout_idx[[i]],i] <- 0
  }
}

png("../figure/experiment/21_camp_dropout.png", height = 800, width = 1500, res = 300, units = "px")
par(mar = rep(0.5,4))
image(.rotate(dropout_mat), breaks = c(-0.5,0.5,1.5), col = c("blue3", "gold"), asp = nrow(dropout_mat)/ncol(dropout_mat),
      axes = F)

for(i in row_idx_lines){
  lines(c(0,1), rep(i, 2), lwd = 3, col = "white")
  lines(c(0,1), rep(i, 2), lwd = 3, lty = 2)
}
graphics.off()

##################

# side-by-side comparison
dat_vec <- as.vector(dat)
dropout_vec <- as.vector(dropout_mat)

x_seq <- seq(min(dat_vec[dat_vec > 0]), max(dat), length.out = 50)
tmp <- diff(x_seq[1:2])/2
x_seq <- x_seq - tmp
x_seq <- c(-tmp, x_seq, x_seq+tmp)
zz <- hist(dat_vec[dropout_vec == 0], breaks = x_seq, plot = F)
zz$density <- log(zz$density + 1)
yy <- hist(dat_vec[dropout_vec == 1], breaks = x_seq, plot = F)
yy$density <- log(yy$density + 1)

png("../figure/experiment/21_camp_dropout_hist.png", height = 1000, width = 2500, res = 300, units = "px")
par(mfrow = c(1,2), mar = c(4,4,3,0.5))
plot(yy, col = "gold", ylim = c(0, max(zz$density)), xlim = c(-.1,2),
     ylab = "Logged density", xlab = "Value", main = "Zoomed in X-axis")
plot(zz, col = rgb(0,0,205/255,0.5), add = T)
legend("topright", c("Dropped out","Kept"),
       bty="n", fill=c("blue3", "gold"));

plot(yy, col = "gold", ylim = c(0, 0.04), xlim = c(0,5),
     ylab = "Logged density", xlab = "Value", main = "Zoomed in Y-axis")
plot(zz, col = rgb(0,0,205/255,0.5), add = T)
graphics.off()

#####################

# plot the weights of zeros across different rows
zero_weight_vec <- rep(NA, ncol(dat))
for(i in 1:length(zero_weight_vec)){
  if(!any(is.na(dropout_res[[i]]))){
    set.seed(10)
    x <- dat[,i]
    x2 <- .jitter_zeros(x)
    val <- min(x2)
    wt <- .calculate_weight(val, dropout_res[[i]])
    zero_weight_vec[i] <- wt[,2]/(wt[,1] + wt[,2])
  }
}

png("../figure/experiment/21_camp_zero_weight.png", height = 1000, width = 1200, res = 300, units = "px")
plot(sort(zero_weight_vec), main = "Probability ratio of 0",
     ylab = "Probability ratio", xlab = "Sorted index of genes",
     pch = 16, col = rgb(0,0,0,0.2))
graphics.off()
length(which(zero_weight_vec >= 0.1))/length(zero_weight_vec)
which(zero_weight_vec >= 0.1)
zero_weight_vec[which(zero_weight_vec >= 0.1)]

png("../figure/experiment/21_camp_hist.png", height = 2000, width = 2500, res = 300, units = "px")
par(mfrow = c(2,2), mar = c(4,3,3,0.5))
idx <- which.max(zero_weight_vec) #3
x <- dat[,idx]
.hist_augment(x, param_list = list(dropout_res[[idx]]), multiplier = 10, lwd = 3, xlab = "Value",
              ylab = "Count", main = paste0("Gene ", idx, ", Ratio at 0: ",
                                            round(zero_weight_vec[idx]*100), "%"))

idx <- 356
x <- dat[,idx]
.hist_augment(x, param_list = list(dropout_res[[idx]]), multiplier = 5, lwd = 3, xlab = "Value",
              ylab = "Count", main = paste0("Gene ", idx, ", Ratio at 0: ",
                                            round(zero_weight_vec[idx]*100), "%"))

idx <- 281
x <- dat[,idx]
.hist_augment(x, param_list = list(dropout_res[[idx]]), multiplier = 0.05, lwd = 3, xlab = "Value",
              ylab = "Count", main = paste0("Gene ", idx, ", Ratio at 0: ",
                                            round(zero_weight_vec[idx]*100), "%"))

idx <- which.min(zero_weight_vec)
x <- dat[,idx]
.hist_augment(x, param_list = list(dropout_res[[idx]]), multiplier = 10, lwd = 3, xlab = "Value",
              ylab = "Count", main = paste0("Gene ", idx, ", Ratio at 0: ",
                                            round(zero_weight_vec[idx]*100), "%"))
graphics.off()

###########################

.find_neighbor <- function(mat, i){
  row_vec <- c(1:nrow(mat))[-i]
  idx_me <- which(mat[i,] == 1)
  vec <- sapply(row_vec, function(x){
    idx_other <- which(mat[x,] == 1)
    length(intersect(idx_me, idx_other))/length(unique(c(idx_me, idx_other)))
  })

  row_vec[order(vec, decreasing = T)[1:25]]
}

neighbor_list <- lapply(1:nrow(dropout_mat), function(x){
  .find_neighbor(dropout_mat, x)
})

.find_neighbor_value <- function(mat, i){
  row_vec <- c(1:nrow(mat))[-i]
  idx_me <- which(mat[i,] == 1)
  vec <- sapply(row_vec, function(x){
    idx_other <- which(mat[x,] == 1)
    length(intersect(idx_me, idx_other))/length(unique(c(idx_me, idx_other)))
  })

  sort(vec, decreasing = T)[25]
}

quantile(sapply(1:nrow(dropout_mat), function(x){
  .find_neighbor_value(dropout_mat, x)
}))

# determine which entries should be zero
zero_mat <- dropout_mat
for(i in 1:nrow(zero_mat)){
  idx <- which(dat[i,] == 0)
  vec <- rep(2, length(idx))
  zz <- colSums(dat[neighbor_list[[i]],idx])
  vec[which(zz == 0)] <- 3
  zero_mat[i, idx] <- vec
}

table(as.numeric(zero_mat))

png("../figure/experiment/21_camp_truezeros.png", height = 800, width = 1500, res = 300, units = "px")
par(mar = rep(0.5,4))
image(.rotate(zero_mat), breaks = c(-0.5,0.5,1.5,2.5,3.5), col = c("blue3", "gold", "blue3", "firebrick1"), asp = nrow(zero_mat)/ncol(zero_mat),
      axes = F)

for(i in row_idx_lines){
  lines(c(0,1), rep(i, 2), lwd = 3, col = "white")
  lines(c(0,1), rep(i, 2), lwd = 3, lty = 2)
}
graphics.off()

# locate the column with most "true" zeros
true_count <- apply(zero_mat, 2, function(x){length(which(x == 3))})
idx <- which.max(true_count)
x <- dat[,idx]
.hist_augment(x, param_list = list(dropout_res[[idx]]), multiplier = 2, lwd = 3, xlab = "Value",
              ylab = "Count")
zero_weight_vec[idx]

# locate the "middle ground"
mid_count <- apply(zero_mat, 2, function(x){length(which(x == 2))/length(which(x == 3))})
idx <- which.min(abs(mid_count - 1))
x <- dat[,idx]
.hist_augment(x, param_list = list(dropout_res[[idx]]), multiplier = 2, lwd = 3, xlab = "Value",
              ylab = "Count")
zero_weight_vec[idx]

# locate the column with most "dropped" zeros
true_count <- apply(zero_mat, 2, function(x){length(which(x == 2))})
idx <- which.max(true_count)
x <- dat[,idx]
.hist_augment(x, param_list = list(dropout_res[[idx]]), multiplier = 2, lwd = 3, xlab = "Value",
              ylab = "Count")
zero_weight_vec[idx]

### see how often it goes outside its cell group
label_vector <- as.numeric(as.factor(sapply(camp$cell.info[,2], function(x){substr(x,start = 1, stop = 1)})))
label_mat <- t(sapply(1:3, function(x){
  idx <- which(label_vector == x)
  neigh_lab <- as.vector(sapply(idx, function(y){
    tmp <- label_vector[neighbor_list[[y]]]
    tmp
  }))
  sapply(1:3, function(k){length(which(neigh_lab == k))})
}))

percentage_mat <- label_mat
for(i in 1:3){
  percentage_mat[i,] <- percentage_mat[i,]/sum(percentage_mat[i,])
}

png("../figure/experiment/21_camp_neighbor.png", height = 1100, width = 1000, res = 300, units = "px")
col_vec <- colorRampPalette(c("white", rgb(0.803, 0.156, 0.211)))(20)
graphics::par(mar = c(1,1,3,0))
image(.rotate(percentage_mat), col = col_vec, asp = T, xlab = "", ylab ="", axes = F,
      main = "Neighboring cell types")
graphics::title(xlab="Neighbor's cell-type", line=0)
graphics::title(ylab="True cell-type", line=0)

for(i in 1:3){
  for(j in 1:3){
    graphics::text((j-1)*.5, 1-(i-1)*0.5,
                   labels = paste0("n = ", label_mat[i,j], "\n(", round(percentage_mat[i,j]*100), "%)"))
  }
}
graphics.off()


#########################

png("../figure/experiment/21_camp_umat.png", height = 800, width = 2000, res = 300, units = "px")
par(mfrow = c(1,3), mar = c(4,4,3,0.5))
# svd with nothing special
res_svd <- svd(dat)
k <- 4
u_mat <- res_svd$u[,1:k] %*% diag(sqrt(res_svd$d[1:k]))
v_mat <- res_svd$u[,1:k] %*% diag(sqrt(res_svd$d[1:k]))

plot(u_mat[,1], u_mat[,2], asp = T, col = label_vector, pch = 16,
     main = "No dropouts", ylab = "U[,2]", xlab = "U[,1]")
lines(rep(0, 2), c(-1e6, 1e6), col = "red", lwd = 2, lty = 2)
lines(c(-1e6, 1e6), rep(0, 2), col = "red", lwd = 2, lty = 2)

# only elements deemed by mixture
dat2 <- dat
dat2[dropout_mat == 0] <- NA
# dat2 <- methods::new("realRatingMatrix", data = recommenderlab::dropNA(dat2))
# funk_svd <- recommenderlab::funkSVD(dat2, k = 4, lambda = 0.01)
# u_mat <- funk_svd$U
# v_mat <- funk_svd$V
softImpute::lambda0(dat2)
set.seed(10)
res1 <- softImpute::softImpute(dat2, rank = 4, lambda = 1)
dim(res1$u)
u_mat <- res1$u %*% diag(sqrt(res1$d))

plot(u_mat[,1], u_mat[,2], xlim = range(as.numeric(u_mat[,1]), 0),
     ylim = range(as.numeric(u_mat[,2]), 0),
     asp = T, col = label_vector, pch = 16,
     main = "Dropouts via mixture model", ylab = "U[,2]", xlab = "U[,1]")
lines(rep(0, 2), c(-1e6, 1e6), col = "red", lwd = 2, lty = 2)
lines(c(-1e6, 1e6), rep(0, 2), col = "red", lwd = 2, lty = 2)

# including elements we now consider true
dat2 <- dat
dat2[dropout_mat == 0] <- NA
dat2[zero_mat == 3] <- 0
# dat2 <- methods::new("realRatingMatrix", data = recommenderlab::dropNA(dat2))
# funk_svd <- recommenderlab::funkSVD(dat2, k = 4, lambda = 0.01)
# u_mat <- funk_svd$U
# v_mat <- funk_svd$V
softImpute::lambda0(dat2)
set.seed(10)
res2 <- softImpute::softImpute(dat2, rank=4, lambda = 1)
u_mat <- res2$u %*% diag(sqrt(res2$d))

plot(u_mat[,1], u_mat[,2], xlim = range(as.numeric(u_mat[,1]), 0),
     ylim = range(as.numeric(u_mat[,2]), 0),
     asp = T, col = label_vector, pch = 16,
     main = "Dropouts via mixture model\nand imputation", ylab = "U[,2]", xlab = "U[,1]")
lines(rep(0, 2), c(-1e6, 1e6), col = "red", lwd = 2, lty = 2)
lines(c(-1e6, 1e6), rep(0, 2), col = "red", lwd = 2, lty = 2)
graphics.off()

# determine the fits
fit_svd <- res_svd$u[,1:k] %*% diag(res_svd$d[1:k]) %*% t(res_svd$v[,1:k])
fit_res1 <- res1$u %*% diag(res1$d) %*% t(res1$v)
fit_res2 <- res2$u %*% diag(res2$d) %*% t(res2$v)

png("../figure/experiment/21_camp_fitted_vals.png", height = 800, width = 2000, res = 300, units = "px")
par(mfrow = c(1,3), mar = c(4,4,3,0.5))
ylim <- range(c(fit_svd, fit_res1, fit_res2))
plot(NA, ylim = ylim, xlim = range(dat), asp = T,
     main = "No dropouts", xlab = "Observed value", ylab = "Predicted value")
lines(rep(0, 2), c(-1e6, 1e6), col = "red", lwd = 1, lty = 1)
lines(c(-1e6, 1e6), rep(0, 2), col = "red", lwd = 1, lty = 1)
lines(c(-1e6, 1e6), c(-1e6, 1e6), col = "red", lwd = 2, lty = 2)
points(as.numeric(dat), as.numeric(fit_svd), pch = 16, col = rgb(0,0,0,0.1))

plot(NA, ylim = ylim, xlim = range(dat), asp = T,
     main = "Dropouts via mixture model", xlab = "Observed value", ylab = "Predicted value")
lines(rep(0, 2), c(-1e6, 1e6), col = "red", lwd = 1, lty = 1)
lines(c(-1e6, 1e6), rep(0, 2), col = "red", lwd = 1, lty = 1)
lines(c(-1e6, 1e6), c(-1e6, 1e6), col = "red", lwd = 2, lty = 2)
points(as.numeric(dat), as.numeric(fit_res1), pch = 16, col = rgb(0,0,0,0.1))

plot(NA, ylim = ylim, xlim = range(dat), asp = T,
     main = "Dropouts via mixture model\nand imputation", xlab = "Observed value", ylab = "Predicted value")
lines(rep(0, 2), c(-1e6, 1e6), col = "red", lwd = 1, lty = 1)
lines(c(-1e6, 1e6), rep(0, 2), col = "red", lwd = 1, lty = 1)
lines(c(-1e6, 1e6), c(-1e6, 1e6), col = "red", lwd = 2, lty = 2)
points(as.numeric(dat), as.numeric(fit_res2), pch = 16, col = rgb(0,0,0,0.1))

graphics.off()

##
png("../figure/experiment/21_camp_fitted_vals_zoom.png", height = 800, width = 2000, res = 300, units = "px")
par(mfrow = c(1,3), mar = c(4,4,3,0.5))
ylim <- range(c(fit_svd, fit_res1, fit_res2))
plot(NA, ylim = c(0, 0.5), xlim = c(0, 0.5), asp = T,
     main = "No dropouts", xlab = "Observed value", ylab = "Predicted value")
lines(rep(0, 2), c(-1e6, 1e6), col = "red", lwd = 1, lty = 1)
lines(c(-1e6, 1e6), rep(0, 2), col = "red", lwd = 1, lty = 1)
lines(c(-1e6, 1e6), c(-1e6, 1e6), col = "red", lwd = 2, lty = 2)
points(as.numeric(dat), as.numeric(fit_svd), pch = 16, col = rgb(0,0,0,0.05))

plot(NA, ylim = c(0, 0.5), xlim = c(0, 0.5), asp = T,
     main = "Dropouts via mixture model", xlab = "Observed value", ylab = "Predicted value")
lines(rep(0, 2), c(-1e6, 1e6), col = "red", lwd = 1, lty = 1)
lines(c(-1e6, 1e6), rep(0, 2), col = "red", lwd = 1, lty = 1)
lines(c(-1e6, 1e6), c(-1e6, 1e6), col = "red", lwd = 2, lty = 2)
points(as.numeric(dat), as.numeric(fit_res1), pch = 16, col = rgb(0,0,0,0.05))

plot(NA, ylim = c(0, 0.5), xlim = c(0, 0.5), asp = T,
     main = "Dropouts via mixture model\nand imputation", xlab = "Observed value", ylab = "Predicted value")
lines(rep(0, 2), c(-1e6, 1e6), col = "red", lwd = 1, lty = 1)
lines(c(-1e6, 1e6), rep(0, 2), col = "red", lwd = 1, lty = 1)
lines(c(-1e6, 1e6), c(-1e6, 1e6), col = "red", lwd = 2, lty = 2)
points(as.numeric(dat), as.numeric(fit_res2), pch = 16, col = rgb(0,0,0,0.05))

graphics.off()

## plot a histogram
idx2 <- which(zero_mat == 2)
idx3 <- which(zero_mat == 3)

png("../figure/experiment/21_camp_hist_zero.png", height = 800, width = 2000, res = 300, units = "px")
par(mfrow = c(1,3), mar = c(4,4,3,0.5))
# x_range <- range(c(fit_svd[c(idx2,idx3)], fit_res1[c(idx2,idx3)], fit_res2[c(idx2,idx3)]))
# x_seq <- seq(x_range[1], x_range[2], length.out = 50)
lower <- -0.05; upper <- 0.25
x_seq <- seq(lower, upper, length.out = 50)
zz2 <- fit_svd[idx2]; zz2 <- zz2[intersect(which(zz2 >= lower), which(zz2 <= upper))]
zz3 <- fit_svd[idx3]; zz3 <- zz3[intersect(which(zz3 >= lower), which(zz3 <= upper))]
hist(zz2, breaks = x_seq, col = "blue3",
     ylab = "Counts", xlab = "Predicted value for 0's",
     main = "No dropouts")
hist(zz3, col = "gold", add = T, breaks = x_seq)
lines(rep(0, 2), c(-1e6,1e6), col = "red", lwd = 2, lty = 2)
legend("topright", c("Dropped zero","Kept zero"),
       bty="n", fill=c("blue3", "gold"))

zz2 <- fit_res1[idx2]; zz2 <- zz2[intersect(which(zz2 >= lower), which(zz2 <= upper))]
zz3 <- fit_res1[idx3]; zz3 <- zz3[intersect(which(zz3 >= lower), which(zz3 <= upper))]
hist(zz2, breaks = x_seq, col = "blue3",
     ylab = "Counts", xlab = "Predicted value for 0's",
     main = "Dropouts via mixture model")
hist(zz3, col = "gold", add = T, breaks = x_seq)
lines(rep(0, 2), c(-1e6,1e6), col = "red", lwd = 2, lty = 2)

zz2 <- fit_res2[idx2]; zz2 <- zz2[intersect(which(zz2 >= lower), which(zz2 <= upper))]
zz3 <- fit_res2[idx3]; zz3 <- zz3[intersect(which(zz3 >= lower), which(zz3 <= upper))]
hist(zz2, breaks = x_seq, col = "blue3",
     ylab = "Counts", xlab = "Predicted value for 0's",
     main = "Dropouts via mixture model\nand imputation",
     ylim = c(0, 5000))
hist(zz3, col = "gold", add = T, breaks = x_seq)
lines(rep(0, 2), c(-1e6,1e6), col = "red", lwd = 2, lty = 2)
graphics.off()

##########



