rm(list=ls())
set.seed(10)
p <- 0.2
dat <- rnbinom(10000, 100, 1-p)
pred_mat <- rep(mean(dat), length(dat))

target_val <- sum((dat - pred_mat)^2)/length(dat)
r_seq <- round(seq(1, 200, length.out = 101))
proposed_val <- sapply(r_seq, function(x){sum((pred_mat + pred_mat^2/x))/length(dat)})
supposed_val1 <- var(dat)
pred_mat + pred_mat^2/100

plot(r_seq, proposed_val, ylim = range(c(target_val, proposed_val)))
lines(c(-1e3, 1e3), rep(target_val, 2), lwd = 2, lty = 2, col = "red")

r_seq[which.min(abs(target_val - proposed_val))]
