% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/initialization.R
\name{.svd_projection}
\alias{.svd_projection}
\title{Do an SVD projection}
\usage{
.svd_projection(mat, k, factors = F, u_alone = F, v_alone = F)
}
\arguments{
\item{mat}{numeric matrix with \code{n} rows and \code{p} columns}

\item{k}{positive integer less than \code{min(c(n,p))}.}

\item{factors}{boolean. If \code{TRUE}, return the factors (i.e., a list containing
\code{u_mat} and \code{v_mat}). If \code{FALSE}, return the rank-\code{k} matrix directly.}

\item{u_alone}{boolean. If \code{TRUE}, then place the singular values on \code{v_mat}}

\item{v_alone}{boolean. If \code{TRUE}, then place the singular values on \code{u_mat}}
}
\value{
a list or numeric matrix, depending on \code{factors}
}
\description{
Uses \code{RSpectra::svds} to compute the \code{k} leading singular vectors, but
sometimes there are numerical instability issues. In case of crashes, the code
then uses the default \code{svd} function.
}
