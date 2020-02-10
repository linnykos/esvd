rm(list-ls())
delta = 0.1
A = matrix(c(1+delta, 0, 0, 1-delta), 2, 2)
H = matrix(c(-delta, delta, delta, delta), 2, 2)
A2 = matrix(c(1,delta,delta,1), 2, 2)

A2 - (A+H)

eigen1 = eigen(A)
eigenH = eigen(H)
eigen2 = eigen(A2)

###############

lhs = t(eigen2$vectors[,2])%*%eigen1$vectors[,1]

rhs_denom = abs(eigen2$values[2] - eigen1$values[1])
rhs_num = t(eigen2$vectors[,2]) %*%H %*% eigen1$vectors[,1]

rhs_num/rhs_denom

##########################3

A = matrix(c(1,0,0,0), 2, 2)
B = matrix(c(1,0,0,0), 2, 2)
AB = A+B

eigen(A)
eigen(B)
eigen(AB)

