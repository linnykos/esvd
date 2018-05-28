load("../../SOUP/data/camp.rda")

dat <- camp$counts

idx <- which(colnames(dat) %in% camp$select.genes)
dat <- dat[,idx]

dat <- t(apply(dat, 1, function(x){x/sum(x)}))
dat <- 1000*log(dat + 1)

dim(dat)

##################

res_svd <- svd(dat)

plot(res_svd$d[1:50])
k <- 4
u_mat <- res_svd$u[,1:k] %*% diag(sqrt(res_svd$d[1:k]))
v_mat <- res_svd$v[,1:k] %*% diag(sqrt(res_svd$d[1:k]))

plot(u_mat[,1], u_mat[,2], asp = T, xlim = range(c(u_mat[,1], 0)), ylim = range(c(u_mat[,2], 0)))
plot(v_mat[,1], v_mat[,2], asp = T, xlim = range(c(v_mat[,1], 0)), ylim = range(c(v_mat[,2], 0)))

u_mat_spherical <- t(apply(u_mat, 1, function(x){x/.l2norm(x)}))
v_mat_spherical <- t(apply(v_mat, 1, function(x){x/.l2norm(x)}))

set.seed(10)
u_clust <- kmeans(u_mat_spherical, centers = 3, iter.max = 100, nstart = 10)
v_clust <- kmeans(v_mat_spherical, centers = 3, iter.max = 100, nstart = 10)

#################

cell_type <- as.numeric(as.factor(camp$cell.info[,2]))
row_idx <- unlist(lapply(1:max(cell_type), function(x){
  which(cell_type == x)
}))
col_idx <- unlist(lapply(1:max(u_clust$cluster), function(x){
  sort(which(v_clust$cluster == x))
}))

dat3 <- dat[row_idx, col_idx]

tmp <- as.numeric(dat3)
tmp <- tmp[tmp!=0]

break_vec <- quantile(tmp, probs = seq(0, 1, length.out = 20))
break_vec <- c(-5, break_vec)

colorRamp_custom <- function(vec1, vec2, length){
  mat <- matrix(0, nrow = length, ncol = 3)
  for(i in 1:3){
    mat[,i] <- seq(vec1[i], vec2[i], length.out = length)
  }

  luminosity_vec <- apply(mat, 1, function(x){
    0.2126*x[1] + 0.7152*x[2] + 0.0722*x[3]
  })

  target_luminosity <- mean(c(luminosity_vec[1], luminosity_vec[length]))

  mat <- t(sapply(1:nrow(mat), function(x){
    factor <- min(c(target_luminosity/luminosity_vec[x], 1/mat[x,]))
    mat[x,] * factor
  }))

  apply(mat, 1, function(x){
    rgb(x[1], x[2], x[3])
  })
}

col_vec <- colorRamp_custom(c(0.584, 0.858, 0.564), c(0.803, 0.156, 0.211), 19)
col_vec <- c("white", col_vec)

.rotate <- function(mat){t(mat)[,nrow(mat):1]}

png("../figure/experiment/6_camp_data_sorted.png", height = 2400, width = 2400, res = 300, units = "px")
image(.rotate(dat3), breaks = break_vec, col = col_vec, asp = nrow(dat3)/ncol(dat3),
      axes = F)

#put lines
row_idx_lines <- sapply(1:(max(cell_type)-1), function(x){
  1-length(which(cell_type <= x))/length(cell_type)
})
col_idx_lines <- sapply(1:(max(v_clust$cluster)-1), function(x){
  length(which(v_clust$cluster <= x))/length(v_clust$cluster)
})

for(i in row_idx_lines){
  lines(c(0,1), rep(i, 2), lwd = 2, lty = 2)
}
for(i in col_idx_lines){
  lines(rep(i, 2), c(0,1), lwd = 2, lty = 2)
}

graphics.off()

######
cell_type <- camp$cell.info[,2]
cell_type_coarse <- as.numeric(as.factor(sapply(cell_type, function(x){substr(x, 1, 1)})))
row_idx <- unlist(lapply(1:max(cell_type_coarse), function(x){
  which(cell_type_coarse == x)
}))

dat4 <- dat[row_idx, col_idx]

png("../figure/experiment/6_camp_data_sorted_3group.png", height = 2400, width = 2400, res = 300, units = "px")
image(.rotate(dat4), breaks = break_vec, col = col_vec, asp = nrow(dat4)/ncol(dat4),
      axes = F)

#put lines
row_idx_lines <- sapply(1:(max(cell_type_coarse)-1), function(x){
  1-length(which(cell_type_coarse <= x))/length(cell_type_coarse)
})
col_idx_lines <- sapply(1:(max(v_clust$cluster)-1), function(x){
  length(which(v_clust$cluster <= x))/length(v_clust$cluster)
})

for(i in row_idx_lines){
  lines(c(0,1), rep(i, 2), lwd = 2, lty = 2)
}
for(i in col_idx_lines){
  lines(rep(i, 2), c(0,1), lwd = 2, lty = 2)
}

graphics.off()

##############################

# crude estimate of dropout
pred_dat <- u_mat %*% t(v_mat)
pred_dat3 <- pred_dat[row_idx, col_idx]
# pred_dat <- res_svd$u[,1:k] %*% diag(res_svd$d[1:k]) %*% t(res_svd$v[,1:k])

png("../figure/experiment/6_camp_data_predicted.png", height = 2400, width = 2400, res = 300, units = "px")
image(.rotate(pred_dat3), breaks = break_vec, col = col_vec, asp = nrow(dat3)/ncol(dat3),
      axes = F)

#put lines
row_idx <- sapply(1:(max(cell_type)-1), function(x){
  1-length(which(cell_type <= x))/length(cell_type)
})
col_idx <- sapply(1:(max(v_clust$cluster)-1), function(x){
  length(which(v_clust$cluster <= x))/length(v_clust$cluster)
})

for(i in row_idx){
  lines(c(0,1), rep(i, 2), lwd = 2, lty = 2)
}
for(i in col_idx){
  lines(rep(i, 2), c(0,1), lwd = 2, lty = 2)
}

graphics.off()

##### now do the calculation

idx <- which(dat == 0)
zero_val <- pred_dat[idx]
hist(zero_val, breaks = 50, col = "gray")

###############################

#let's see how good naive k-means clustering does

set.seed(10)
k_clust <- kmeans(dat, centers = 3, nstart = 100, iter.max = 100)
table(k_clust$cluster, cell_type_coarse)

##############################


