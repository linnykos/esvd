rm(list = ls())
load("../experiment/lpsolve_error/example.RData")
res <- lpSolve::lp("min", objective_in, const_mat, const_dir, const_rhs)

# no feasible solution found, however:

zz <- const_mat %*% c(current_vec, 0, 0) # <- this is the allowable solution that lpSolve misses

k <- nrow(const_mat)/2
all(zz[1:k] == zz[(k+1):(2*k)])

all(zz[1:k] <= const_rhs[1:k]) &
  all(zz[(k+1):(2*k)] >= const_rhs[(k+1):(2*k)])
all(const_dir[1:k] == "<=")
all(const_dir[(k+1):(2*k)] == ">=")
