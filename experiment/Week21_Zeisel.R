rm(list=ls())
load("../../SOUP/data/zeisel.rda")

dat <- zeisel$counts

dat <- t(apply(dat, 1, function(x){x/sum(x)}))
dat <- 100*log10(dat + 1)

idx <- which(colnames(dat) %in% zeisel$select.genes)
dat <- dat[,idx]
dim(dat)

table(zeisel$cell.info$cell.type)
idx <- grep("oligodendrocytes", zeisel$cell.info$cell.type)
dat_subset <- dat[idx,]
dim(dat_subset)

##########################################################

dropout_res <- vector("list", length = ncol(dat_subset))
for(i in 1:ncol(dat_subset)){
  set.seed(10)
  if(i %% floor(ncol(dat_subset)/10) == 0) cat('*')
  x <- dat_subset[,i]
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
weight_mat <- matrix(NA, ncol = ncol(dat_subset), nrow = nrow(dat_subset))
for(i in 1:ncol(weight_mat)){
  if(any(is.na(dropout_res[[i]]))){
    weight_mat[,i] <- 0
  } else {
    set.seed(10)
    x <- dat_subset[,i]
    x2 <- .jitter_zeros(x)
    wt <- .calculate_weight(x2, dropout_res[[i]])
    weight_mat[,i] <- wt[,2]/(wt[,1] + wt[,2])
  }
}

.rotate <- function(mat){t(mat)[,nrow(mat):1]}

col_vec <- colorRampPalette(c("blue3", "gold"))(20)
png("../figure/experiment/21_zeisel_weight.png", height = 800, width = 2400, res = 300, units = "px")
par(mar = rep(0.5,4))
image(.rotate(weight_mat), col = col_vec, asp = nrow(weight_mat)/ncol(weight_mat),
      axes = F)
graphics.off()

dropout_idx <- vector("list", length = ncol(dat_subset))
for(i in 1:ncol(dat_subset)){
  if(!any(is.na(dropout_res[[i]]))){
    set.seed(10)
    x <- dat_subset[,i]
    x2 <- .jitter_zeros(x)
    dropout_idx[[i]] <- which(.compute_dropout(dropout_res[[i]], x2) == 1)
  }
}

#1 means keep, 0 means drop
dropout_mat <- matrix(1, ncol = ncol(dat_subset), nrow = nrow(dat_subset))
for(i in 1:ncol(dat_subset)){
  if(length(dropout_idx[[i]]) == 0){
    dropout_mat[,i] <- 0
  } else {
    dropout_mat[dropout_idx[[i]],i] <- 0
  }
}

png("../figure/experiment/21_zeisel_dropout.png", height = 800, width = 2400, res = 300, units = "px")
par(mar = rep(0.5,4))
image(.rotate(dropout_mat), breaks = c(-0.5,0.5,1.5), col = c("blue3", "gold"), asp = nrow(dropout_mat)/ncol(dropout_mat),
      axes = F)
graphics.off()

########################################

# side-by-side comparison
dat_vec <- as.vector(dat_subset)
dropout_vec <- as.vector(dropout_mat)

x_seq <- seq(min(dat_vec[dat_vec > 0]), max(dat_vec), length.out = 50)
tmp <- diff(x_seq[1:2])/2
x_seq <- x_seq - tmp
x_seq <- c(-tmp, x_seq, x_seq+tmp)
zz <- hist(dat_vec[dropout_vec == 0], breaks = x_seq, plot = F)
zz$density <- log(zz$density + 1)
yy <- hist(dat_vec[dropout_vec == 1], breaks = x_seq, plot = F)
yy$density <- log(yy$density + 1)

png("../figure/experiment/21_zeisel_dropout_hist.png", height = 1000, width = 2500, res = 300, units = "px")
par(mfrow = c(1,2), mar = c(4,4,3,0.5))
plot(zz, col = "gold", ylim = c(0, max(yy$density)), xlim = c(-.1,2),
     ylab = "Logged density", xlab = "Value", main = "Zoomed in X-axis")
plot(yy, col = rgb(0,0,205/255,0.5), add = T)
legend("topright", c("Dropped out","Kept"),
       bty="n", fill=c("blue3", "gold"));

plot(zz, col = "gold", ylim = c(0, 0.04), xlim = c(0,5),
     ylab = "Logged density", xlab = "Value", main = "Zoomed in Y-axis")
plot(yy, col = rgb(0,0,205/255,0.5), add = T)
graphics.off()

##################################

# plot the weights of zeros across different rows
zero_weight_vec <- rep(NA, ncol(dat_subset))
for(i in 1:length(zero_weight_vec)){
  if(!any(is.na(dropout_res[[i]]))){
    x <- dat_subset[,i]
    x2 <- .jitter_zeros(x)
    val <- min(x2)
    wt <- .calculate_weight(val, dropout_res[[i]])
    zero_weight_vec[i] <- wt[,2]/(wt[,1] + wt[,2])
  }
}
plot(sort(zero_weight_vec))
length(which(zero_weight_vec >= 0.1))/length(zero_weight_vec)

idx <- which.max(zero_weight_vec) #1195
x <- dat_subset[,idx]
.hist_augment(x, param_list = list(dropout_res[[idx]]), multiplier = 2, lwd = 3, xlab = "Value",
              ylab = "Count")

order(zero_weight_vec, decreasing = T)[1:100]
idx <- 510
x <- dat_subset[,idx]
.hist_augment(x, param_list = list(dropout_res[[idx]]), multiplier = 2, lwd = 3, xlab = "Value",
              ylab = "Count")
zero_weight_vec[510]

idx <- 1330
x <- dat_subset[,idx]
.hist_augment(x, param_list = list(dropout_res[[idx]]), multiplier = 2, lwd = 3, xlab = "Value",
              ylab = "Count")
zero_weight_vec[1330]

#########################

# let's try this idea of seeing if the nearest neighbors also have a zero
## to define a nearest neighbor, we see the percentage overlap
.find_neighbor <- function(mat, i){
  row_vec <- c(1:nrow(mat))[-i]
  idx_me <- which(mat[i,] == 1)
  vec <- sapply(row_vec, function(x){
    idx_other <- which(mat[x,] == 1)
    length(intersect(idx_me, idx_other))/length(unique(c(idx_me, idx_other)))
  })

  row_vec[order(vec, decreasing = T)[1:50]]
}

neighbor_list <- lapply(1:nrow(dropout_mat), function(x){
  print(x)
  .find_neighbor(dropout_mat, x)
})

# determine which entries should be zero
zero_mat <- dropout_mat
for(i in 1:nrow(zero_mat)){
  print(i)
  idx <- which(dat_subset[i,] == 0)
  vec <- rep(2, length(idx))
  zz <- colSums(dat_subset[neighbor_list[[i]],idx])
  vec[which(zz == 0)] <- 3
  zero_mat[i, idx] <- vec
}

table(as.numeric(zero_mat))

png("../figure/experiment/21_zeisel_truezeros.png", height = 800, width = 2400, res = 300, units = "px")
par(mar = rep(0.5,4))
image(.rotate(zero_mat), breaks = c(-0.5,0.5,1.5,2.5,3.5), col = c("blue3", "gold", "blue3", "firebrick1"), asp = nrow(zero_mat)/ncol(zero_mat),
      axes = F)
graphics.off()

# locate the column with most "true" zeros
true_count <- apply(zero_mat, 2, function(x){length(which(x == 3))})
discard_idx <- which(sapply(dropout_res, function(x){any(is.na(x))}))
idx <- c(1:ncol(zero_mat))[-discard_idx][which.max(true_count[-discard_idx])]
x <- dat_subset[,idx]
.hist_augment(x, param_list = list(dropout_res[[idx]]), multiplier = 2, lwd = 3, xlab = "Value",
              ylab = "Count")
zero_weight_vec[idx]

# locate the "middle ground"
mid_count <- apply(zero_mat, 2, function(x){length(which(x == 2))/length(which(x == 3))})
idx <- which.min(abs(mid_count - 1))
x <- dat_subset[,idx]
.hist_augment(x, param_list = list(dropout_res[[idx]]), multiplier = 2, lwd = 3, xlab = "Value",
              ylab = "Count")
zero_weight_vec[idx]

#########

# TRY THIS WITH CAMP
# see if the eigenvectors change dramatically
