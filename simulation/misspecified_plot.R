rm(list=ls())
load("../results/misspecified_simulation.RData")

col_vec <- c(rgb(205,40,54,maxColorValue=255), #red
             rgb(180,200,255,maxColorValue=255), #blue
             rgb(100,100,200,maxColorValue=255), #purple
             rgb(149,219,144,maxColorValue=255)) #green

for(i in 1:length(res)){
  plot(res[[i]][[1]]$u_mat[,1], res[[i]][[1]]$u_mat[,2], pch = 16,
       asp = T, col = col_vec[rep(1:4, each = paramMat[1,"n"])],
       main = i)
}
