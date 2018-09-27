rm(list=ls())
source("../experiment/Week23_simulation_generator.R")
library(singlecell)

col_vec <- c(rgb(205,40,54,maxColorValue=255),
             rgb(180,200,255,maxColorValue=255),
             rgb(100,100,200,maxColorValue=255),
             rgb(149,219,144,maxColorValue=255))

set.seed(10)
simulation <- .data_generator(total = 150, distr_func = function(x){stats::rnorm(1, x, x/2)})
dat <- simulation$dat
length(which(dat == 0))/prod(dim(dat))
.plot_singlecell(dat)

par(mfrow = c(1,3))
zz <- apply(dat, 1, function(x){length(which(x != 0))/length(x)})
plot(sort(zz))
zz <- apply(dat, 2, function(x){length(which(x != 0))/length(x)})
plot(sort(zz))
plot(sort(dat[dat != 0]))


###########

# try the naive thing
svd_res <- svd(dat)
k <- 2
u_mat <- svd_res$u %*% diag(sqrt(svd_res$d))
library(slingshot)
slingshot_res <- slingshot::slingshot(data = u_mat, clusterLabels = rep(1:simulation$h, each = simulation$n_each))

plot(u_mat[,1], u_mat[,2],
     col = col_vec[rep(1:simulation$h, each = simulation$n_each)], asp = T,
     pch = 16, xlab = "Estimated dim. 1", ylab = "Estimated dim. 2", main = "Cell estimated vectors")
lines(slingshot_res)

##############

# let's do something a bit more sensible?
dropout_mat <- singlecell:::.dropout(dat)
table(dropout_mat)

zero_mat <- singlecell:::.find_true_zeros(dropout_mat)
c(table(zero_mat), length(which(is.na(zero_mat))))

dat2 <- dat
dat2[which(is.na(zero_mat))] <- NA
res <- singlecell:::.fit_gaussian_factorization(dat2)
