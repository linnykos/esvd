\name{LLE}
\alias{LLE}
\title{
	Locally Linear Embedding
}
\description{
	Computes the Locally Linear Embedding as introduced in 2000 by Roweis, Saul and Lawrence. 
}
\usage{
LLE(data, dim=2, k)
}
\arguments{
  \item{data}{
	N x D matrix (N samples, D features)
}
  \item{dim}{
	dimension of the target space
}
  \item{k}{
	number of neighbours
}
}
\details{

	Locally Linear Embedding (LLE) preserves local properties of the data by
	representing each sample in the data by a linear combination of
	its k nearest neighbours with each neighbour weighted
	independently. LLE finally chooses the low-dimensional
	representation that best preserves the weights in the target
	space. \cr
	This R version is based on the Matlab implementation by Sam Roweis.
}
\value{
	It returns a N x dim matrix (N samples, dim features) with the reduced input data
}
\references{
	Roweis, Sam T. and Saul, Lawrence K., "Nonlinear Dimensionality Reduction by Locally Linear Embedding",2000;
}
\author{
	Christoph Bartenhagen
}