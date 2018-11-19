library(descend)
data(zeisel)

result <- runDescend(zeisel$count.matrix.small,
                     scaling.consts = zeisel$library.size, n.cores = 1)

hvg <- findHVG(result)

which(hvg$score.mat[,2]>2)
