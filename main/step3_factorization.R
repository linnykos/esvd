load("../results/step1_imputing.RData")

max_val <- 1000
scalar_vec <- c(0.1, 0.25, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 8, 10, 100)
res_list <- vector("list", length(scalar_vec))

mean_val <- mean(dat_impute)
dat_impute <- t(apply(dat_impute, 1, function(x){x/sum(x)}))
dat_impute <- dat_impute * mean_val/mean(dat_impute)

for(i in 1:length(scalar_vec)){
  init <- singlecell:::.initialization(dat_impute, family = "gaussian", scalar = scalar_vec[i],
                                       k = 5, max_val = max_val)
  res_list[[i]] <- singlecell:::.fit_factorization(dat_impute, u_mat = init$u_mat, v_mat = init$v_mat,
                                                   family = "gaussian",  reparameterize = T,
                                                   max_iter = 25, max_val = max_val,
                                                   scalar = scalar_vec[i],
                                                   return_path = F, cores = 15,
                                                   verbose = T)
  save.image("../results/step3_factorization_tmp.RData")
}

save.image("../results/step3_factorization.RData")
