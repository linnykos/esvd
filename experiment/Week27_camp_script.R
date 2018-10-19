rm(list=ls())
load("../../SOUPR/data/camp.rda")

dat <- camp$counts

idx <- which(colnames(dat) %in% camp$select.genes)
dat <- dat[,idx]

dat <- t(apply(dat, 1, function(x){x/sum(x)}))
dat <- 10*log(dat + 1)

dim(dat)
dat <- dat[order(camp$cell.info[,2]),]
camp$cell.info <- camp$cell.info[order(camp$cell.info[,2]),]

plot(sort(as.numeric(dat[dat > 0])))
yy <- sort(as.numeric(dat[dat > 0]))

# zz <- dat[dat > 0]
# max_val <- -1/mean(zz[zz < quantile(zz, probs = 0.2)])
max_val <- -1/yy[7500]
7500/prod(dim(dat))
max_iter <- 50

###############

dropout_mat <- singlecell:::.dropout(dat)
zero_mat <- singlecell:::.find_true_zeros(dropout_mat, num_neighbors = 50)
length(which(zero_mat == 0))/prod(dim(zero_mat))

set.seed(10)
dat_impute <- singlecell:::.scImpute(dat, which(is.na(zero_mat)), Kcluster = 7,
                                     verbose = T, weight = 0.5)

init_impute <- singlecell:::.initialization(dat_impute, family = "exponential",
                                            max_val = max_val)

set.seed(10)
res_withimpute <- singlecell:::.fit_factorization(dat_impute, init_impute$u_mat, init_impute$v_mat,
                                                  verbose = T, family = "exponential",
                                                  max_iter = max_iter, tol = NA,
                                                  cores = 15, max_val = max_val)

save.image("Week27_camp.RData")



