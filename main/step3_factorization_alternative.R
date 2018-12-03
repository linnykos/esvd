load("../results/step1_imputing.RData")

max_val <- 1000
scalar_vec <- c(0.1, 0.25, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 8, 10, 100)
res_list <- vector("list", length(scalar_vec))

mean_val <- mean(dat_impute)
dat_impute <- t(apply(dat_impute, 1, function(x){x/sum(x)}))
dat_impute <- dat_impute * mean_val/mean(dat_impute)


## generate some missing values

set.seed(10)
idx <- rbind(do.call(rbind, (lapply(1:n, function(x){
  cbind(x, sample(1:d, 4))
}))), do.call(rbind, (lapply(1:d, function(x){
  cbind(sample(1:n, 4), d)
}))))

dat_impute_NA <- dat_impute
for(i in 1:nrow(idx)){
  dat_impute_NA[idx[i,1], idx[i,2]] <- NA
}
idx <- which(is.na(dat_impute_NA))

res_naive <- softImpute::softImpute(dat_impute_NA, rank.max = 5, lambda = 30)
pred_naive <- res_naive$u %*% diag(res_naive$d) %*% t(res_naive$v)

# plot(pred_mat[idx], dat_impute[idx], pch = 16, asp = T, col = rgb(0,0,0,0.2))
# lines(c(-1e10,1e10), c(-1e10, 1e10), col = "red", lwd = 2)

##################

max_val <- 1000
scalar_vec <- c(0.1, 0.25, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 8, 10, 100)
res_list <- vector("list", length(scalar_vec))

mean_val <- mean(dat_impute)
dat_impute <- t(apply(dat_impute, 1, function(x){x/sum(x)}))
dat_impute <- dat_impute * mean_val/mean(dat_impute)

for(i in 1:length(scalar_vec)){
  init <- singlecell:::.initialization(dat_impute_NA, family = "gaussian", scalar = scalar_vec[i],
                                       k = 5, max_val = max_val)
  res_list[[i]] <- singlecell:::.fit_factorization(dat_impute_NA, u_mat = init$u_mat, v_mat = init$v_mat,
                                                   family = "gaussian",  reparameterize = T,
                                                   max_iter = 25, max_val = max_val,
                                                   scalar = scalar_vec[i],
                                                   return_path = F, cores = 15,
                                                   verbose = T)
  save.image("../results/step3_factorization_NA_tmp.RData")
}

save.image("../results/step3_factorization_NA.RData")
