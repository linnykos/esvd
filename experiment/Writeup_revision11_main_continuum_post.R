rm(list=ls())
load("../results/Writeup_revision11_continuum.RData")

par(mfrow = c(1,2))
nrow_vec <- c(nrow(dat1), nrow(dat2))
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

  plot(NA, ylim = range(obj_vec_list[[k]]), xlim = c(0, nrow_vec[k]),
       main = paste0(nrow_vec[k]))
  for(i in 1:length(start_vec_list[[k]])){
    lines(x = c(start_vec_list[[k]][i], end_vec_list[[k]][i]), y = rep(obj_vec_list[[k]][i], 2), lwd = 2)
  }
}

########################

zz <- intersect(which(obj_vec_list[[1]]>2), which(midpoint_vec_list[[1]] >= 2000))
head(zz)

j <- zz[3]
par(mfrow = c(1,2))
vec1 <- dat1[,j]
col_vec <- rep("black", length(vec1)); col_vec[segmentation_res[[j]]$cut_1$i:segmentation_res[[j]]$cut_1$j] <- "red"
vec2 <- dat2[,j]
plot(vec1, col = col_vec, pch = 16, cex = 0.5)
lines(rep(max_common_idx, 2), c(-1e5, 1e5), col = "red", lwd = 2)
points(segmentation_res[[j]]$vec1_smooth, col = "white", cex = 0.75, pch = 16)
points(segmentation_res[[j]]$vec1_smooth, col = "red", cex = 0.5, pch = 16)
plot(vec2, pch = 16, cex = 0.5)
lines(rep(max_common_idx, 2), c(-1e5, 1e5), col = "red", lwd = 2)
points(segmentation_res[[j]]$vec2_smooth, col = "white", cex = 0.75, pch = 16)
points(segmentation_res[[j]]$vec2_smooth, col = "red", cex = 0.5, pch = 16)

#######

zz <- intersect(which(obj_vec_list[[2]]>2), which(midpoint_vec_list[[1]] >= 2000))
head(zz)

j <- zz[2]
par(mfrow = c(1,2))
vec1 <- dat1[,j]
vec2 <- dat2[,j]
col_vec <- rep("black", length(vec2)); col_vec[segmentation_res[[j]]$cut_2$i:segmentation_res[[j]]$cut_2$j] <- "red"
plot(vec1, pch = 16, cex = 0.5)
lines(rep(max_common_idx, 2), c(-1e5, 1e5), col = "red", lwd = 2)
points(segmentation_res[[j]]$vec1_smooth, col = "white", cex = 0.75, pch = 16)
points(segmentation_res[[j]]$vec1_smooth, col = "red", cex = 0.5, pch = 16)
plot(vec2, col = col_vec, pch = 16, cex = 0.5)
lines(rep(max_common_idx, 2), c(-1e5, 1e5), col = "red", lwd = 2)
points(segmentation_res[[j]]$vec2_smooth, col = "white", cex = 0.75, pch = 16)
points(segmentation_res[[j]]$vec2_smooth, col = "red", cex = 0.5, pch = 16)


