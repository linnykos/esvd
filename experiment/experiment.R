rm(list=ls())
load("../results/factorization_results_others.RData")

col_func <- function(alpha){
  c( rgb(86/255, 180/255, 233/255, alpha), #skyblue
     rgb(240/255, 228/255, 66/255, alpha), #yellow
     rgb(0/255, 158/255, 115/255, alpha), #bluish green
     rgb(230/255, 159/255, 0/255,alpha)) #orange
}
col <- col_func(1)

cluster_labels <- rep(1:4, each = paramMat[1,"n_each"])

type_vec <- c("gaussian", "curved_gaussian", "zinb", "pcmf")
k <- 1
idx <- 1
i <- 1
curves <- singlecell::slingshot(res[[k]][[idx]][[i]], cluster_labels,
                                starting_cluster = 1, verbose = F,
                                use_initialization = T)

####################

