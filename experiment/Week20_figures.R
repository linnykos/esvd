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

###############

png("../figure/experiment/20_zeisel_illustrative.png", height = 1200, width = 2200, res = 300, units = "px")
par(mfrow = c(2,2), mar = c(4,4,0.5,0.5))
x <- dat_subset[,1302]
.hist_augment(x, xlab = "Value", ylab = "Count", main = "")

x <- dat_subset[,1043]
.hist_augment(x, xlab = "Value", ylab = "Count", main = "")

x <- dat_subset[,340]
.hist_augment(x, xlab = "Value", ylab = "Count", main = "")

x <- dat_subset[,1249]
.hist_augment(x, xlab = "Value", ylab = "Count", main = "")
graphics.off()

############################

rm(list=ls())
load("../../SOUP/data/camp.rda")
dat <- camp$counts

dat <- t(apply(dat, 1, function(x){x/sum(x)}))
dat <- log10(dat + 1)

idx <- which(colnames(dat) %in% camp$select.genes)
dat <- dat[,idx]
dim(dat)

zz <- svd(dat)
k <- 5
dat_pred <- zz$u[,1:k] %*% diag(zz$d[1:k]) %*% t(zz$v[,1:k])

png("../figure/experiment/20_camp_predict.png", height = 1500, width = 1500, res = 300, units = "px")
plot(as.numeric(dat), as.numeric(dat_pred), pch = 16, asp = T,
     xlab = "Observed value", ylab = "Predicted value",
     main = "", col = rgb(0,0,0,0.1), xlim = c(0, 0.0075),
     ylim = c(0, 0.0075))
lines(c(-1e6, 1e6), c(-1e6, 1e6), col = "red", lwd = 1, lty = 2)
lines(c(-1e6, 1e6), rep(0,2), col = "red", lwd = 1)
lines(rep(0,2), c(-1e6, 1e6), col = "red", lwd = 1)
graphics.off()
