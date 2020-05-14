rm(list=ls())
load("../results/tmp_continuum2.RData")

# NOTE TO SELF: I bungled the code, but we're applying circular_list to each column of the dat_impute here

j <- 1

zz <- order(apply(pred_mat, 2, max), decreasing = T)
zz[1:10]
# zz = order(apply(dat_impute, 2, function(x){quantile(x, probs = 0.75)}), decreasing = T)
# zz[1:20]

for(j in zz[1:10]){
  par(mfrow = c(1,2))
  vec <- pred_mat[idx_cell[order(order_vec)], j]
  vec2 <- dat_impute[idx_cell[order(order_vec)], j]
  col_vec <- rep(rgb(0,0,0,0.1), length(vec))
  col_vec[circular_list[[j]]$i : circular_list[[j]]$j] <- "red"
  plot(vec, ylim = range(c(vec, vec2)), col = col_vec, pch = 16,
       main = paste0("Predicted gene expression\n(Gene ", colnames(dat_impute)[j], ")"))
  plot(vec2, ylim = range(c(vec, vec2)), col = col_vec, pch = 16,
       main = paste0("Observed gene expression\n(Gene ", colnames(dat_impute)[j], ")"))
}

############

# now extract the midpoints
start_vec <- sapply(1:length(circular_list), function(i){circular_list[[i]]$i})
end_vec <- sapply(1:length(circular_list), function(i){circular_list[[i]]$j})
midpoint_vec <- sapply(1:length(start_vec), function(i){(start_vec[i] + end_vec[i])/2})
length_vec <- end_vec - start_vec
obj_vec <- sapply(1:length(circular_list), function(i){circular_list[[i]]$obj_val})
gene_vec <- 1:ncol(dat_impute)

plot(NA, ylim = range(log(pmax(obj_vec,0)+1)), xlim = c(0, length(idx_cell)))
for(i in 1:length(start_vec)){
  lines(x = c(start_vec[i], end_vec[i]), y = rep(log(max(obj_vec[i],0)+1), 2), lwd = 2)
}

idx <- which(length_vec >= length(idx_cell)/2)

start_vec <- start_vec[-idx]
end_vec <- end_vec[-idx]
midpoint_vec <- midpoint_vec[-idx]
length_vec <- length_vec[-idx]
obj_vec <- obj_vec[-idx]
gene_vec <- gene_vec[-idx]

plot(NA, ylim = range(log(pmax(obj_vec,0)+1)), xlim = c(0, length(idx_cell)))
for(i in 1:length(start_vec)){
  lines(x = c(start_vec[i], end_vec[i]), y = rep(log(max(obj_vec[i],0)+1), 2), lwd = 2)
}

# select highly expressed genes
selected_idx <- which(log(pmax(obj_vec,0)+1) >= 1)
# order these genes
selected_idx <- selected_idx[order(midpoint_vec[selected_idx], decreasing = F)]

perc <- .1
j <- gene_vec[selected_idx[round(perc * length(selected_idx))]]
par(mfrow = c(1,2))
vec <- pred_mat[idx_cell[order(order_vec)], j]
vec2 <- dat_impute[idx_cell[order(order_vec)], j]
col_vec <- rep(rgb(0,0,0,0.1), length(vec))
col_vec[circular_list[[j]]$i : circular_list[[j]]$j] <- "red"
plot(vec, ylim = range(c(vec, vec2)), col = col_vec, pch = 16,
     main = paste0("Predicted gene expression\n(Gene ", colnames(dat_impute)[j], ")"))
plot(vec2, ylim = range(c(vec, vec2)), col = col_vec, pch = 16,
     main = paste0("Observed gene expression\n(Gene ", colnames(dat_impute)[j], ")"))
