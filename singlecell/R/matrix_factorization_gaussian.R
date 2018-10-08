.evaluate_objective.gaussian <- function(dat, u_mat, v_mat, ...){
  stopifnot(is.matrix(dat), nrow(dat) == nrow(u_mat), ncol(dat) == nrow(v_mat),
            ncol(u_mat) == ncol(v_mat))


}
