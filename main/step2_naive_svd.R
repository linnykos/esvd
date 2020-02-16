set.seed(10)
load(paste0("../results/step1_imputing", suffix, ".RData"))

n <- nrow(dat_impute); d <- ncol(dat_impute)
set.seed(10)
missing_idx <- rbind(do.call(rbind, (lapply(1:n, function(x){
  cbind(x, sample(1:d, 4))
}))), do.call(rbind, (lapply(1:d, function(x){
  cbind(sample(1:n, 4), d)
}))))

dat_impute_NA <- dat_impute
for(i in 1:nrow(missing_idx)){
  dat_impute_NA[missing_idx[i,1], missing_idx[i,2]] <- NA
}
missing_idx <- which(is.na(dat_impute_NA))

# k <- 4
# lambda0_val <- softImpute::lambda0(dat_impute_NA)
# res_naive <- softImpute::softImpute(dat_impute_NA, rank.max = k, lambda = min(30, lambda0_val/100))
# pred_naive <- res_naive$u %*% diag(res_naive$d) %*% t(res_naive$v)

rm(list = c("n", "d", "lambda0_val", "res_hvg"))
print(paste0(Sys.time(), ": Finished naive SVD"))
save.image(paste0("../results/step2_naive_svd", suffix, ".RData"))
