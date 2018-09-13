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

set.seed(10)
png("../figure/experiment/20_zeisel.png", height = 2000, width = 2500, res = 300, units = "px")
par(mfrow = c(2,2), mar = c(4,3,3,0.5))
x <- dat_subset[,17]
zz <- .em_mixture(x)
x2 <- .jitter_zeros(x)
count <- 100*length(which(.compute_dropout(zz, x2) == 1))/length(x)
.hist_augment(x, param_list = list(zz), multiplier = 10, lwd = 3, xlab = "Value",
              ylab = "Count", main = paste0("Gene 17, (", round(count, 1), "% dropped)"))
legend("topright", c("T. Gaussian distribution","Gamma distribution"),
       bty="n", fill=c(rgb(0.584, 0.858, 0.564), rgb(0.803, 0.156, 0.211)));
# zz <- .em_mixture(x, mixture = "exponential.tgaussian")
# .hist_augment(x, param_list = list(zz), multiplier = 0.02, lwd = 3)

x <- dat_subset[,340]
zz <- .em_mixture(x)
x2 <- .jitter_zeros(x)
count <- 100*length(which(.compute_dropout(zz, x2) == 1))/length(x)
.hist_augment(x, param_list = list(zz), multiplier = 100, lwd = 3, xlab = "Value",
              ylab = "Count", main = paste0("Gene 340, (", round(count, 1), "% dropped)"))
# zz <- .em_mixture(x, mixture = "exponential.tgaussian")
# .hist_augment(x, param_list = list(zz), multiplier = 1, lwd = 3)

x <- dat_subset[,1249]
zz <- .em_mixture(x)
count <- 100*length(which(.compute_dropout(zz, x) == 1))/length(x)
.hist_augment(x, param_list = list(zz), multiplier = 100, lwd = 3, xlab = "Value",
              ylab = "Count", main = paste0("Gene 1249, (", round(count, 1), "% dropped)"))

x <- dat_subset[,1900]
zz <- .em_mixture(x)
x2 <- .jitter_zeros(x)
count <- 100*length(which(.compute_dropout(zz, x2) == 1))/length(x)
.hist_augment(x, param_list = list(zz), multiplier = 25, lwd = 3, xlab = "Value",
              ylab = "Count", main = paste0("Gene 1900, (", round(count, 1), "% dropped)"))
# zz <- .em_mixture(x, mixture = "exponential.tgaussian")
# .hist_augment(x, param_list = list(zz), multiplier = 0.1, lwd = 3)
graphics.off()

##################################################

dropout_idx <- vector("list", length = ncol(dat_subset))
for(i in 1:ncol(dat_subset)){
  set.seed(10)
  if(i %% floor(ncol(dat_subset)/10) == 0) cat('*')
  x <- dat_subset[,i]
  if(all(x == 0)) {
    dropout_idx[[i]] <- NA
    next()
  }
  if(length(x[x > 0]) <= 10) {
    dropout_idx[[i]] <- NA
    next()
  }

  dropout_idx[[i]] <- tryCatch({
    res <- .em_mixture(x)
    x2 <- .jitter_zeros(x)
    which(.compute_dropout(res, x2) == 1)
  }, error = function(e){
    print(paste("Error in ", i))
    return(NA)
  })
}

cutoff_vec <- sapply(1:ncol(dat_subset), function(i){
  x <- dat_subset[,i]
  if(all(is.na(x))){
    NA
  } else if (length(dropout_idx[[i]]) == nrow(dat_subset)){
    NA
  } else {
    idx1 <- which(x > 0)
    if(length(idx1) == 0) return(NA)
    idx2 <- c(1:nrow(dat_subset))[-dropout_idx[[i]]]
    if(length(intersect(idx1, idx2)) == 0) return(min(x[idx1]))
    min(x[intersect(idx1, idx2)])
  }
})

min(cutoff_vec, na.rm = T)
zz <- cbind(cutoff_vec, 1:ncol(dat_subset), apply(dat_subset, 2, function(x){length(x[x>0])}))
zz[order(zz[,1], decreasing = F)[1:100],]


png("../figure/experiment/20_zeisel_2.png", height = 1000, width = 2500, res = 300, units = "px")
par(mfrow = c(1,2), mar = c(4,3,3,0.5))
set.seed(10)
x <- dat_subset[,566]
zz <- .em_mixture(x)
x2 <- .jitter_zeros(x)
count <- 100*length(which(.compute_dropout(zz, x2) == 1))/length(x)
.hist_augment(x, param_list = list(zz), multiplier = 10, lwd = 3, xlab = "Value",
              ylab = "Count", main = paste0("Gene ,566, (", round(count, 1), "% dropped)"))
legend("topright", c("T. Gaussian distribution","Gamma distribution"),
       bty="n", fill=c(rgb(0.584, 0.858, 0.564), rgb(0.803, 0.156, 0.211)));

set.seed(10)
x <- dat_subset[,1302]
zz <- .em_mixture(x)
x2 <- .jitter_zeros(x)
count <- 100*length(which(.compute_dropout(zz, x2) == 1))/length(x)
.hist_augment(x, param_list = list(zz), multiplier = 10, lwd = 3, xlab = "Value",
              ylab = "Count", main = paste0("Gene 1302, (", round(count, 1), "% dropped)"))
graphics.off()



count_vec <- sapply(dropout_idx, function(x){
  if(all(is.na(x))){
    NA
  } else {
    length(x)
  }
})

quantile(count_vec, na.rm = T)
order(count_vec, decreasing = F)[1:2]

png("../figure/experiment/20_zeisel_3.png", height = 1000, width = 2500, res = 300, units = "px")
par(mfrow = c(1,2), mar = c(4,3,3,0.5))
set.seed(10)
x <- dat_subset[,1043]
zz <- .em_mixture(x)
x2 <- .jitter_zeros(x)
count <- 100*length(which(.compute_dropout(zz, x2) == 1))/length(x)
.hist_augment(x, param_list = list(zz), multiplier = 10, lwd = 3, xlab = "Value",
              ylab = "Count", main = paste0("Gene 1043, (", round(count, 1), "% dropped)"))
legend("topright", c("T. Gaussian distribution","Gamma distribution"),
       bty="n", fill=c(rgb(0.584, 0.858, 0.564), rgb(0.803, 0.156, 0.211)));

set.seed(10)
x <- dat_subset[,1175]
zz <- .em_mixture(x)
x2 <- .jitter_zeros(x)
count <- 100*length(which(.compute_dropout(zz, x2) == 1))/length(x)
.hist_augment(x, param_list = list(zz), multiplier = 10, lwd = 3, xlab = "Value",
              ylab = "Count", main = paste0("Gene 1175, (", round(count, 1), "% dropped)"))
graphics.off()

#####################################

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

tmp <- as.numeric(dat_subset)
tmp <- tmp[tmp!=0]

break_vec <- quantile(tmp, probs = seq(0, 1, length.out = 20))
break_vec <- c(-5, break_vec)

png("../figure/experiment/20_zeisel_data.png", height = 800, width = 2400, res = 300, units = "px")
par(mar = rep(0.5,4))
image(.rotate(dat_subset), breaks = break_vec, col = col_vec, asp = nrow(dat_subset)/ncol(dat_subset),
      axes = F)
graphics.off()

dropout_mat <- matrix(0, ncol = ncol(dat_subset), nrow = nrow(dat_subset))
for(i in 1:ncol(dat_subset)){
  dropout_mat[dropout_idx[[i]],i] <- 1
}

png("../figure/experiment/20_zeisel_dropout.png", height = 800, width = 2400, res = 300, units = "px")
par(mar = rep(0.5,4))
image(.rotate(dropout_mat), breaks = c(-0.5,0.5,1.5), col = c("gold", "blue3"), asp = nrow(dropout_mat)/ncol(dropout_mat),
      axes = F)
graphics.off()

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

png("../figure/experiment/20_zeisel_dropout_hist.png", height = 1000, width = 2500, res = 300, units = "px")
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
