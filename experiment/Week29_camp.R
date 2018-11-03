rm(list=ls())
load("../../SOUPR/data/camp.rda")

dat <- camp$counts

idx <- which(colnames(dat) %in% camp$select.genes)
dat <- dat[,idx]
dim(dat)

set.seed(0)
# res <- SAVER::saver(dat)
# dim(res$estimate)
# save.image("../experiment/Week29_camp_saver_impute.RData")

load("../experiment/Week29_camp_saver_impute.RData")

cell_type <- camp$cell.info[,2]
cell_type_coarse <- as.numeric(as.factor(sapply(cell_type, function(x){substr(x, 1, 1)})))
cell_type_numeric <- as.numeric(as.factor(cell_type))

col_vec <- c(rgb(205,40,54,maxColorValue=255), #red
             rgb(180,200,255,maxColorValue=255), #purple
             rgb(100,100,200,maxColorValue=255), #blue
             rgb(149,219,144,maxColorValue=255)) #green

###################################################


svd_res <- svd(res$estimate)
k <- 8
u_mat_naive <- svd_res$u[,1:k] %*% diag(sqrt(svd_res$d[1:k]))
v_mat_naive <- svd_res$v[,1:k] %*% diag(sqrt(svd_res$d[1:k]))

png("../figure/experiment/29_camp_native.png", height = 1400, width = 2400, res = 300, units = "px")
par(mfrow = c(1,2))
plot(-u_mat_naive[,1], u_mat_naive[,2],
     col = col_vec[cell_type_coarse],
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
