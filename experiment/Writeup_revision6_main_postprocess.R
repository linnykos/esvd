rm(list=ls())
load("../results/step5_clustering.RData")
zz = prcomp(dat_impute, center = T, scale. = F)
vec <- zz$rotation[,1]
vec[vec <= 0] <- 0
vec <- (vec - min(vec))/(max(vec) - min(vec))
vec <- vec/sum(vec)
xx <- dat_impute %*% vec

color_func <- function(alpha = 0.2){
  c(rgb(240/255, 228/255, 66/255, alpha), #yellow
    rgb(86/255, 180/255, 233/255, alpha), #skyblue
    rgb(0/255, 158/255, 115/255, alpha), #bluish green
    rgb(0/255, 114/255, 178/255,alpha), #blue
    rgb(230/255, 159/255, 0/255,alpha), #orange
    rgb(210/255, 198/255, 36/255, alpha)) #darker yellow
}
cell_type_vec <- as.character(marques$cell.info$cell.type[cell_idx])
cell_type_vec <- as.factor(cell_type_vec)
cluster_labels <- as.numeric(cell_type_vec)

custom_cluster_group_list <- list(13, 12, 1, c(10,11), c(2,3), c(4:7), c(8,9))

col_vec <- color_func(1)[c(5, rep(3,2), rep(1,6), rep(2,2),  rep(5,2))]
col_name <- c("orange", rep("bluish green", 2), rep("yellow", 6), rep("skyblue", 2), rep("orange", 2))
col_info <- data.frame(name = levels(cell_type_vec),
                       idx = sort(unique(cluster_labels)),
                       level = sapply(1:13, function(y){which(sapply(cluster_group_list, function(x){y %in% x}))}),
                       col_name = col_name,
                       col_code = col_vec)
col_info

tmp_df <- data.frame(val = xx, type = sapply(1:13, function(y){which(sapply(cluster_group_list, function(x){y %in% x}))})[cluster_labels])

vioplot::vioplot(tmp_df$val[tmp_df$type == 1],
                 tmp_df$val[tmp_df$type == 2],
                 tmp_df$val[tmp_df$type == 3],
                 tmp_df$val[tmp_df$type == 4],
                 tmp_df$val[tmp_df$type == 5],
                 tmp_df$val[tmp_df$type == 6],
                 col = sapply(1:6, function(x){unique(col_info$col_code[which(col_info$level == x)])}),
                 pchMed = 21,
                 colMed = "black", colMed2 = "white",
                 xlab = "", names = rep("", 6))

#######################################


