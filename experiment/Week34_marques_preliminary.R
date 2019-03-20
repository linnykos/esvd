rm(list=ls())
load("../results/step3_scalar_heuristic.RData")
plot(quality_vec) # uh oh... not looking great...

# svd (huh... doesn't look that bad)
plot(pred_naive[missing_idx], dat_impute[missing_idx], asp = T, pch = 16,
     col = rgb(0,0,0,0.2))
lines(c(-1e10,1e10), c(-1e10,1e10), col = "red", lwd = 2)

idx <- which.min(quality_vec)
for(i in 1:length(res_list)){
  pred_mat <- 1/(res_list[[i]]$u_mat %*% t(res_list[[i]]$v_mat))
  plot(pred_mat[missing_idx], dat_impute[missing_idx], asp = T, main = i,pch = 16,
       col = rgb(0,0,0,0.2))
  lines(c(-1e10,1e10), c(-1e10,1e10), col = "red", lwd = 2)
}



###########################

rm(list=ls())
load("../results/step4_factorization.RData")
cell_type_vec <- as.character(marques$cell.info$cell.type[cell_idx])
# keep only the first two letters
cell_type_vec <- as.character(sapply(cell_type_vec, function(x){substr(x,1,2)}))
col_idx <- c(rgb(0,0,0,0.3), rgb(1,0,0,0.3), rgb(0,1,0,0.3), rgb(0,0,1,0.3),
             rgb(0,0,0,0.3), rgb(0,0,0,0.3))[as.numeric(as.factor(cell_type_vec))]

num_cell <- length(unique(cell_type_vec))
i1 <- 1; i2 <- 2
xlim <- range(res_our$u_mat[,i1])
ylim <- range(res_our$u_mat[,i2])
par(mfrow = c(2,3))
for(i in 1:6){
  idx <- which(as.numeric(as.factor(cell_type_vec)) == i)
  plot(res_our$u_mat[idx,i1], res_our$u_mat[idx,i2], asp = T, pch = 16,
       col = col_idx[idx], xlim = xlim, ylim = ylim)

  lines(c(-100,100), rep(0,2), lwd = 2, lty = 2)
  lines(rep(0,2), c(-100,100), lwd = 2, lty = 2)
}

#all in the same plot
plot(res_our$u_mat[,1], res_our$u_mat[,2], asp = T, pch = 16,
     col = col_idx)


