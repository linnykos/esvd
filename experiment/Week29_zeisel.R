rm(list=ls())
load("../../SOUPR/data/zeisel.rda")

dat <- zeisel$counts

idx <- which(colnames(dat) %in% zeisel$select.genes)
dat <- dat[,idx]
dim(dat)

table(zeisel$cell.info$cell.type)
idx <- grep("oligodendrocytes", zeisel$cell.info$cell.type)
dat <- dat[idx,]
dim(dat)

# screen out the genes with almost all 0's
zz <- apply(dat, 2, function(x){length(which(x != 0))/length(x)})
idx <- which(zz >= 0.05)
dat <- dat[,idx]
dim(dat)

# .plot_singlecell(dat)

#set.seed(0)
#res <- SAVER::saver(dat)
#save.image("../experiment/Week29_zeisel_saver_impute.RData")
load("../experiment/Week29_zeisel_saver_impute.RData")


svd_res <- svd(res$estimate)
k <- 5
u_mat_naive <- svd_res$u[,1:k] %*% diag(sqrt(svd_res$d[1:k]))
v_mat_naive <- svd_res$v[,1:k] %*% diag(sqrt(svd_res$d[1:k]))

png("../figure/experiment/29_zeisel_native.png", height = 1400, width = 2400, res = 300, units = "px")
par(mfrow = c(1,2))
plot(-u_mat_naive[,1], u_mat_naive[,2],
     xlim = range(c(-u_mat_naive[,1], 0)),
     ylim = range(c(u_mat_naive[,2], 0)),
     pch = 16, xlab = "Latent dim. 1", ylab = "Latent dim. 2", main = "Cell latent vectors\n(Naive fit, imputed mean)")
lines(c(-1e6, 1e6), rep(0, 2), col = "red", lwd = 2, lty = 2)
lines( rep(0, 2), c(-1e6, 1e6), col = "red", lwd = 2, lty = 2)

plot(-v_mat_naive[,1], v_mat_naive[,2],
     xlim = range(c(-v_mat_naive[,1], 0)),
     ylim = range(c(v_mat_naive[,2], 0)),
     pch = 16, xlab = "Latent dim. 1", ylab = "Latent dim. 2", main = "Gene latent vectors\n(Naive fit, imputed mean)")
lines(c(-1e6, 1e6), rep(0, 2), col = "red", lwd = 2, lty = 2)
lines( rep(0, 2), c(-1e6, 1e6), col = "red", lwd = 2, lty = 2)
graphics.off()

####################################

# log scale
range(res$estimate)
log_dat <- log(res$estimate)
svd_res <- svd(log_dat)
k <- 5
u_mat_naive <- svd_res$u[,1:k] %*% diag(sqrt(svd_res$d[1:k]))
v_mat_naive <- svd_res$v[,1:k] %*% diag(sqrt(svd_res$d[1:k]))

png("../figure/experiment/29_zeisel_log.png", height = 1400, width = 2400, res = 300, units = "px")
par(mfrow = c(1,2))
plot(-u_mat_naive[,1], u_mat_naive[,2],
     xlim = range(c(-u_mat_naive[,1], 0)),
     ylim = range(c(u_mat_naive[,2], 0)),
     pch = 16, xlab = "Latent dim. 1", ylab = "Latent dim. 2", main = "Cell latent vectors\n(Naive fit, imputed mean)")
lines(c(-1e6, 1e6), rep(0, 2), col = "red", lwd = 2, lty = 2)
lines( rep(0, 2), c(-1e6, 1e6), col = "red", lwd = 2, lty = 2)

plot(-v_mat_naive[,1], v_mat_naive[,2],
     xlim = range(c(-v_mat_naive[,1], 0)),
     ylim = range(c(v_mat_naive[,2], 0)),
     pch = 16, xlab = "Latent dim. 1", ylab = "Latent dim. 2", main = "Gene latent vectors\n(Naive fit, imputed mean)")
lines(c(-1e6, 1e6), rep(0, 2), col = "red", lwd = 2, lty = 2)
lines( rep(0, 2), c(-1e6, 1e6), col = "red", lwd = 2, lty = 2)
graphics.off()

#################


####################

# dat <- t(apply(dat, 1, function(x){x/sum(x)}))
# dat <- 10*log(dat + 1)
# dropout_res <- lapply(1:ncol(dat), function(i){
#   .em_mixture(dat[,i])
# })
# dropout_mat <- singlecell:::.dropout(dat)
# zero_mat <- singlecell:::.find_true_zeros(dropout_mat, num_neighbors = 50)
# set.seed(10)
# dat_impute <- singlecell:::.scImpute(dat, which(is.na(zero_mat)), Kcluster = 7,
#                                      verbose = T, weight = 0.25)

