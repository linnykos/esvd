rm(list=ls())
load("../results/submission_round2b_revision/step7_additional_analyses.RData")

cell_majortype <- as.factor(sapply(cell_type_vec, function(x){
  substr(as.character(x), 0, 3)
}))
tmp1 <- data.frame(cell_majortype = cell_majortype,
                   cell_subtype = cell_type_vec, coord = esvd_embedding$u_mat)

write.csv(tmp1, file = "../results/embedding_coordinates.csv", row.names = F)

curve1 <- esvd_curves_prepared$Curve1
curve2 <- esvd_curves_prepared$Curve2

nrow(curve1)
nrow(curve2)

idx <- c(1,which(sapply(1:(nrow(curve1)-1), function(i){
  sum(abs(curve1[i,] - curve1[i+1,])) >= 1e-6
})), nrow(curve1))
curve1 <- curve1[idx,]

idx <- c(1,which(sapply(1:(nrow(curve2)-1), function(i){
  sum(abs(curve2[i,] - curve2[i+1,])) >= 1e-6
})), nrow(curve2))
curve2 <- curve2[idx,]

nrow(curve1)
nrow(curve2)

curve1 <- data.frame(coord = curve1)
curve2 <- data.frame(coord = curve2)

write.csv(curve1, file = "../results/curve1_coordinates.csv", row.names = F)
write.csv(curve2, file = "../results/curve2_coordinates.csv", row.names = F)
