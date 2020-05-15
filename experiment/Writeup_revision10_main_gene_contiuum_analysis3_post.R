rm(list=ls())
load("../results/Writeup_revision10_main_gene_continuum_analysis3.RData")

par(mfrow = c(1,2))
nrow_vec <- c(nrow(dat_ordered1), nrow(dat_ordered2))
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

  plot(NA, ylim = range(log(pmax(obj_vec_list[[k]],0)+1)), xlim = c(0, nrow_vec[k]))
  for(i in 1:length(start_vec)){
    lines(x = c(start_vec_list[[k]][i], end_vec_list[[k]][i]), y = rep(log(max(obj_vec_list[[k]][i],0)+1), 2), lwd = 2)
  }
}

gene_idx <- intersect(which(log(pmax(obj_vec_list[[1]],0)+1) >= 1),
                      which(midpoint_vec_list[[1]] <= 2200))
zz <- gene_idx[order(obj_vec_list[[1]][gene_idx], decreasing = T)]
head(zz)

j <- 52
par(mfrow = c(1,2))
for(k in 1:2){
  col_vec <- rep("black", nrow(dat_ordered1))
  col_vec[segmentation_res[[j]][[k]]$i : segmentation_res[[j]][[k]]$j] <- "red"
  if(k == 1){
    plot(dat_ordered1[,j], pch = 16, col = col_vec)
    lines(rep(2200, 2), c(-1e5,1e5), col = "red", lwd = 2, lty = 2)
  } else {
    plot(dat_ordered2[,j], pch = 16, col = col_vec)
    lines(rep(2200, 2), c(-1e5,1e5), col = "red", lwd = 2, lty = 2)
  }
}

#######################################

# highly expressed in first trajectory
gene_idx <- intersect(which(log(pmax(obj_vec_list[[1]],0)+1) >= 1),
                      which(midpoint_vec_list[[1]] >= 2500))
zz <- gene_idx[order(obj_vec_list[[1]][gene_idx], decreasing = T)]
head(zz)

j <- 553
par(mfrow = c(1,2))
for(k in 1:2){
  col_vec <- rep("black", nrow(dat_ordered1))
  col_vec[segmentation_res[[j]][[k]]$i : segmentation_res[[j]][[k]]$j] <- "red"
  if(k == 1){
    plot(dat_ordered1[,j], pch = 16, col = col_vec)
    lines(rep(2200, 2), c(-1e5,1e5), col = "red", lwd = 2, lty = 2)
  } else {
    plot(dat_ordered2[,j], pch = 16)
    lines(rep(2200, 2), c(-1e5,1e5), col = "red", lwd = 2, lty = 2)
  }
}

# #### the "consensus" region is the first 2200 indices
# dat_all <- rbind(dat_ordered1[round(seq(1, 2200, length.out = 600)),], dat_ordered1[idx_trajectory1,], dat_ordered2[idx_trajectory2,])
#
#
# j <- 453
# plot(dat_all[,j])


