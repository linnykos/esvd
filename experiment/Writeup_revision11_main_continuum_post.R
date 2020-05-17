rm(list=ls())
load("../results/Writeup_revision11_continuum.RData")

par(mfrow = c(1,2))
nrow_vec <- max_common_idx+c(length(idx_trajectory1), length(idx_trajectory2))
start_vec_list <- vector("list", 2)
end_vec_list <- vector("list", 2)
midpoint_vec_list <- vector("list", 2)
obj_vec_list <- vector("list", 2)

for(k in 1:2){
  start_vec_list[[k]] <- sapply(1:length(segmentation_res), function(i){segmentation_res[[i]][[k]]$i})
  end_vec_list[[k]] <- sapply(1:length(segmentation_res), function(i){segmentation_res[[i]][[k]]$j})
  midpoint_vec_list[[k]] <- sapply(1:length(start_vec_list[[k]]),
                                   function(i){(start_vec_list[[k]][i] + end_vec_list[[k]][i])/2})
  obj_vec_list[[k]] <- sapply(1:length(segmentation_res), function(i){segmentation_res[[i]][[k]]$obj_val})

  plot(NA, ylim = range(log(pmax(obj_vec_list[[k]],0)+1)), xlim = c(0, nrow_vec[k]),
       main = paste0(nrow_vec[k]))
  for(i in 1:length(start_vec_list[[k]])){
    lines(x = c(start_vec_list[[k]][i], end_vec_list[[k]][i]), y = rep(log(max(obj_vec_list[[k]][i],0)+1), 2), lwd = 2)
  }
}

########################

zz <- intersect(which(log(pmax(obj_vec_list[[1]],0)+1)>1.2), which(midpoint_vec_list[[1]] >= 2200))
head(zz)

j <- 184
par(mfrow = c(1,2))
vec1 <- c(dat_ordered1[1:max_common_idx,j], dat_ordered1[idx_trajectory1,j])
col_vec <- rep("black", length(vec1)); col_vec[segmentation_res[[j]]$cut_1$i:segmentation_res[[j]]$cut_1$j] <- "red"
vec2 <- c(dat_ordered1[1:max_common_idx,j], dat_ordered2[idx_trajectory2,j])
plot(vec1, col = col_vec, pch = 16)
lines(rep(max_common_idx, 2), c(-1e5, 1e5), col = "red", lwd = 2)
plot(vec2, pch = 16)
lines(rep(max_common_idx, 2), c(-1e5, 1e5), col = "red", lwd = 2)
