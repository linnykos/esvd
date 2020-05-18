rm(list=ls())
load("../results/step3_baron_factorization_spca_descend.RData")

par(mfrow = c(2,3))
dim1 <- 1; dim2 <- 3
for(i in 1:length(preprocessing_list)){
  n <- nrow(preprocessing_list[[i]]$dat_impute)
  p <- ncol(preprocessing_list[[i]]$dat_impute)
  plot((n/p)^(1/4)*svd_embedding_list[[i]]$u[,dim1]*sqrt(svd_embedding_list[[i]]$d[dim1]),
       (n/p)^(1/4)*svd_embedding_list[[i]]$u[,dim2]*sqrt(svd_embedding_list[[i]]$d[dim2]), asp = T,
       pch = 16, col = as.numeric(as.factor(preprocessing_list[[i]]$label_vec)))
}

for(i in 1:length(preprocessing_list)){
  plot(esvd_embedding_list[[i]]$u[,1], esvd_embedding_list[[i]]$u[,3], asp = T,
       pch = 16, col = as.numeric(as.factor(preprocessing_list[[i]]$label_vec)))
}

# i had this purity idea.... also, it seems like we should be selected the negative binomials...

dat_i <- 1
plot(esvd_missing_list_list[[dat_i]][[8]][[1]]$u_mat[,1],
     esvd_missing_list_list[[dat_i]][[8]][[1]]$u_mat[,2], asp = T,
     pch = 16, col = as.numeric(as.factor(preprocessing_list[[dat_i]]$label_vec)))

# plot(svd_embedding_list[[dat_i]]$u[,1],
#      svd_embedding_list[[dat_i]]$u[,2], asp = T,
#      pch = 16, col = as.numeric(as.factor(preprocessing_list[[dat_i]]$label_vec)))


i <- 6
u_mat <- esvd_embedding_list[[i]]$u_mat
# u_mat <- svd_embedding_list[[i]]$u[,c(1:3)] %*% diag(sqrt(svd_embedding_list[[i]]$d[c(1:3)]))
rgl::plot3d(u_mat[,1], u_mat[,2], u_mat[,3], asp = T, col = as.numeric(preprocessing_list[[i]]$label_vec),
            main = "Dataset 4")
