load("../results/step3_factorization_logged.RData")

k <- 3
u_mat <- res$u_mat[,1:k]
cluster_labels <- dbscan(u_mat, neighbor_count = 10, upper_cutoff = 14,
                     size_cutoff = 19)
# c(length(is.na(cluster_labels)), table(cluster_labels))

upscale_vec <- c(0.5,1,1, 10,10,1, 10,10)
curves <- slingshot(u_mat, cluster_labels, starting_cluster = 1, b = 1, shrink = 1,
                    upscale_vec = upscale_vec)

rm(list = c("u_mat"))
save.image("../results/step4_clustering.RData")
