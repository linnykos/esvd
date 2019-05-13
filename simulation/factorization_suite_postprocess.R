rm(list=ls())
load("../results/factorization_results.RData")

res_mat <- matrix(NA, 6, trials)
for(i in 1:trials){
  for(k in 1:6){
    res_mat[k,i] <- mean(sapply(1:2, function(j){
      idx <- which(res[[1]][[i]]$curves_truth$curves[[j]]$lambda_long != 0)

      # find which curve is most suitable
      tmp <- sapply(1:length(res[[1]][[i]][[2*k]]$curves), function(l){
        cor(res[[1]][[i]]$curves_truth$curves[[j]]$lambda_long[idx],
            res[[1]][[i]][[2*k]]$curves[[l]]$lambda_long[idx], method = "kendall")
      })

      max(tmp)
    }))
  }
}

png("../figure/simulation/factorization_density.png",
    height = 1200, width = 2000, res = 300, units = "px")
label_vec <- c("SVD", "ICA", "t-SNE", "Our", "ZINB-WaVE", "pCMF")
par(mfrow = c(2,3), mar = c(4, 4, 4, 0.5))
for(i in 1:6){
  hist(res_mat[i,], xlim = c(0, 1), breaks = seq(0, 1, length.out = 21), col = "gray",
       xlab = "Kendall's tau correlation", ylab = "Count",
       main = paste0(label_vec[i], ": (", round(median(res_mat[i,]), 2), ")"))
  lines(rep(median(res_mat[i,]), 2), c(0, 100), col = "red", lwd = 2, lty = 2)
}
graphics.off()

############################


col_func <- function(alpha){
  c( rgb(86/255, 180/255, 233/255, alpha), #skyblue
     rgb(240/255, 228/255, 66/255, alpha), #yellow
     rgb(0/255, 158/255, 115/255, alpha), #bluish green
     rgb(230/255, 159/255, 0/255,alpha)) #orange
}
col <- col_func(1)

png("../figure/simulation/factorization_example.png",
    height = 1500, width = 2000, res = 300, units = "px")
label_vec <- c("SVD", "ICA", "t-SNE", "Our", "ZINB-WaVE", "pCMF")
idx <- 50
par(mfrow = c(2,3), mar = c(4, 4, 4, 0.5))
for(i in 1:6){
  plot(res[[1]][[idx]][[(i-1)*2+1]][,1], res[[1]][[idx]][[(i-1)*2+1]][,2],
       asp = T, pch = 16, col = col[rep(1:4, each = paramMat[1,"n_each"])],
       xlab = "Latent dim. 1", ylab = "Latent dim. 2",
       main = paste0(label_vec[i], ": (", round(res_mat[i, idx], 2), ")"))
  for(j in 1:length(res[[1]][[idx]][[(i-1)*2+2]]$curves)){
    ord <- res[[1]][[idx]][[(i-1)*2+2]]$curves[[j]]$ord
    lines(res[[1]][[idx]][[(i-1)*2+2]]$curves[[j]]$s[ord,1],
          res[[1]][[idx]][[(i-1)*2+2]]$curves[[j]]$s[ord,2], lwd = 3, col = "white")
    lines(res[[1]][[idx]][[(i-1)*2+2]]$curves[[j]]$s[ord,1],
          res[[1]][[idx]][[(i-1)*2+2]]$curves[[j]]$s[ord,2], lwd = 2)}
}
graphics.off()

###################


col_func <- function(alpha){
  c( rgb(86/255, 180/255, 233/255, alpha), #skyblue
     rgb(240/255, 228/255, 66/255, alpha), #yellow
     rgb(0/255, 158/255, 115/255, alpha), #bluish green
     rgb(230/255, 159/255, 0/255,alpha)) #orange
}
col <- col_func(1)

paramMat <- cbind(50, 120, 0.05, 150, 2, 2, -2000)
colnames(paramMat) <- c("n_each", "d_each", "sigma", "total", "k", "scalar", "max_val")
vec <- paramMat[1,]
cell_pop <- matrix(c(4,10, 25,100, 60,80, 25,100,
                     40,10, 60,80, 60,80, 100, 25)/10,
                   nrow = 4, ncol = 4, byrow = T)
gene_pop <- matrix(c(20,90, 25,100,
                     90,20, 100,25)/100, nrow = 2, ncol = 4, byrow = T)
n_each <- vec["n_each"]
d_each <- vec["d_each"]
sigma <- vec["sigma"]
total <- vec["total"]

h <- nrow(cell_pop)
cell_mat <- do.call(rbind, lapply(1:h, function(x){
  pos <- stats::runif(n_each)
  cbind(pos*cell_pop[x,1] + (1-pos)*cell_pop[x,3] + stats::rnorm(n_each, sd = sigma),
        pos*cell_pop[x,2] + (1-pos)*cell_pop[x,4] + stats::rnorm(n_each, sd = sigma))
}))
n <- nrow(cell_mat)
k <- ncol(cell_mat)

# construct the gene information
g <- nrow(gene_pop)
gene_mat <- do.call(rbind, lapply(1:g, function(x){
  pos <- stats::runif(d_each)
  cbind(pos*gene_pop[x,1] + (1-pos)*gene_pop[x,3] + stats::rnorm(d_each, sd = sigma),
        pos*gene_pop[x,2] + (1-pos)*gene_pop[x,4] + stats::rnorm(d_each, sd = sigma))
}))
d <- nrow(gene_mat)

# form observations
gram_mat <- cell_mat %*% t(gene_mat) #natural parameter
svd_res <- svd(gram_mat)
cell_mat <- svd_res$u[,1:k] %*% diag(sqrt(svd_res$d[1:k]))
gene_mat <- svd_res$v[,1:k] %*% diag(sqrt(svd_res$d[1:k]))

tmp <- singlecell:::.reparameterize(cell_mat, gene_mat)
cell_mat <- tmp$u_mat; gene_mat <- tmp$v_mat

png("../figure/simulation/factorization_truth.png",
    height = 1200, width = 1200, res = 300, units = "px")
par(mar = c(4, 4, 4, 0.5))
plot(cell_mat[,1], cell_mat[,2],
     asp = T, pch = 16, col = col[rep(1:4, each = paramMat[1,"n_each"])],
     xlab = "Latent dim. 1", ylab = "Latent dim. 2",
     main = "True embedding")
for(j in 1:length(res[[1]][[1]]$curves_truth$curves)){
  ord <- res[[1]][[1]]$curves_truth$curves[[j]]$ord
  lines(res[[1]][[1]]$curves_truth$curves[[j]]$s[ord,1],
        res[[1]][[1]]$curves_truth$curves[[j]]$s[ord,2], lwd = 6, col = "white")
  lines(res[[1]][[1]]$curves_truth$curves[[j]]$s[ord,1],
        res[[1]][[1]]$curves_truth$curves[[j]]$s[ord,2], lwd = 4)}
graphics.off()


