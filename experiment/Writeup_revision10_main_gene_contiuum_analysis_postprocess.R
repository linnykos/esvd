rm(list=ls())
load("../results/tmp_continuum.RData")

zz <- order(apply(pred_mat, 2, max), decreasing = T)
# zz = order(apply(dat_impute, 2, function(x){quantile(x, probs = 0.75)}), decreasing = T)

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

########################

# now extract the midpoints
midpoint_vec <- sapply(1:length(circular_list), function(i){
  (circular_list[[i]]$i + circular_list[[i]]$j)/2
})
start_vec <- sapply(1:length(circular_list), function(i){circular_list[[i]]$i})
end_vec <- sapply(1:length(circular_list), function(i){circular_list[[i]]$j})
plot(sort(midpoint_vec))
obj_vec <- sapply(1:length(circular_list), function(i){circular_list[[i]]$obj_val})
# plot(midpoint_vec, log(pmax(obj_vec,0)+1))
plot(NA, ylim = range(log(pmax(obj_vec,0)+1)), xlim = c(0, length(idx_cell)))
for(i in 1:length(start_vec)){
  lines(x = c(start_vec[i], end_vec[i]), y = rep(log(max(obj_vec[i],0)+1), 2), lwd = 2)
}
