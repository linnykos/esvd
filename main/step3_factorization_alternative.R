load("../results/step1_imputing.RData")

max_val <- 1000
scalar_vec <- c(0.1, 0.25, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 8, 10, 100)
res_list <- vector("list", length(scalar_vec))



## generate some missing values

# plot(pred_mat[idx], dat_impute[idx], pch = 16, asp = T, col = rgb(0,0,0,0.2))
# lines(c(-1e10,1e10), c(-1e10, 1e10), col = "red", lwd = 2)

##################


save.image("../results/step3_factorization_NA.RData")
