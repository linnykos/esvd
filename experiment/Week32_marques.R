rm(list=ls())
load("../results/step3_scalar_heuristic.RData")
plot(scalar_vec, quality_vec) #monkas

rm(list=ls())
load("../results/step4_factorization.RData")
cell_type_vec <- as.character(marques$cell.info$cell.type[cell_idx])
# keep only the first two letters
cell_type_vec <- as.character(sapply(cell_type_vec, function(x){substr(x,1,2)}))
col_idx <- as.numeric(as.factor(cell_type_vec))

num_cell <- length(unique(cell_type_vec))
i1 <- 2; i2 <- 3
xlim <- range(res_our$u_mat[,i1])
ylim <- range(res_our$u_mat[,i2])
par(mfrow = c(2,3))
for(i in 1:6){
  idx <- which(col_idx == i)
  plot(res_our$u_mat[idx,i1], res_our$u_mat[idx,i2], asp = T, pch = 16,
       col = col_idx[idx], xlim = xlim, ylim = ylim)

  lines(c(-100,100), rep(0,2), lwd = 2, lty = 2)
  lines(rep(0,2), c(-100,100), lwd = 2, lty = 2)
}

