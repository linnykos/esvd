% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/matrix_factorization.R
\name{.check_rank}
\alias{.check_rank}
\title{Check the rank of the inner product between the rows of two matrices}
\usage{
.check_rank(u_mat, v_mat)
}
\arguments{
\item{u_mat}{a matrix}

\item{v_mat}{a matrix with the same number of columns as \code{u_mat}}
}
\value{
a list containing \code{u_mat} and \code{v_mat}
}
\description{
If the rank of \code{u_mat \%*\% t(v_mat)} is not equal to
\code{ncol(u_mat)} (which is equal to \code{ncol(v_mat)}), then
drop the last few columns of \code{u_mat} and \code{v_mat}.
This function assumes the source of lack-of-rank is from these dropped columns.
}
