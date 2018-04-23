rm(list=ls())

statements <- vector("list")

load("/raid6/Kevin/singlecell_results/mice/Zeisel_expr.RData")

dat <- as.matrix(dat)

statements[[1]] <- dim(dat)
statements[[2]] <- table(as.vector(dat))

load("/raid6/Kevin/singlecell_results/mice/svd_tmp.RData")

statements[[3]] <- res$d
statements[[4]] <- res$u[,1:3]
statements[[5]] <- res$v[,1:3]

res2 <- res
res2$d <- res$d[1:6]
res2$v <- res$v[,1:6]
res2$u <- res$u[,1:6]

dat_approx <- res2$u %*% diag(res2$d) %*% t(res2$v)
residual <- dat - dat_approx

set.seed(10)
size <- 1e6
idx <- sample(1:prod(dim(dat)), size)
df <- data.frame(mean = as.vector(dat_approx)[idx], resid = as.vector(residual)[idx])

statements[[6]] <- df

load("/raid6/Kevin/singlecell_results/mice/isomap_tmp.RData")

statements[[7]] <- isomap$dim2

#############

load("/raid6/Kevin/singlecell_results/mice/Zeisel_expr.RData")

dat <- as.matrix(dat)
logdat <- log(dat+1)

load("/raid6/Kevin/singlecell_results/mice/svd_tmp_log.RData")

statements[[8]] <- res$d
statements[[9]] <- res$u[,1:2]
statements[[10]] <- res$v[,1:2]

res2 <- res
res2$d <- res$d[1:2]
res2$v <- res$v[,1:2]
res2$u <- res$u[,1:2]

logdat_approx <- res2$u %*% diag(res2$d) %*% t(res2$v)
residual <- logdat - logdat_approx

set.seed(10)
idx <- sample(1:prod(dim(dat)), size)
df <- data.frame(mean = as.vector(logdat_approx)[idx], resid = as.vector(residual)[idx])

statements[[11]] <- df

###############


load("/raid6/Kevin/singlecell_results/mice/Zeisel_expr.RData")

dat <- as.matrix(dat)
dat <- t(apply(dat, 1, function(x){x/sum(x)}))
logdat <- log(dat+1)

load("/raid6/Kevin/singlecell_results/mice/svd_tmp_log_renormalized.RData")

statements[[12]] <- res$d
statements[[13]] <- res$u[,1:2]
statements[[14]] <- res$v[,1:2]

res2 <- res
res2$d <- res$d[1:2]
res2$v <- res$v[,1:2]
res2$u <- res$u[,1:2]

logdat_approx <- res2$u %*% diag(res2$d) %*% t(res2$v)
residual <- logdat - logdat_approx

set.seed(10)
idx <- sample(1:prod(dim(dat)), size)
df <- data.frame(mean = as.vector(logdat_approx)[idx], resid = as.vector(residual)[idx])

statements[[15]] <- df

##############

load("/raid6/Kevin/singlecell_results/mice/Zeisel_expr.RData")

dat <- as.matrix(dat)
sum_vec <- apply(dat, 1, sum)
statements[[16]] <- sum_vec


save(statements, file = "statements_tmp.RData")


