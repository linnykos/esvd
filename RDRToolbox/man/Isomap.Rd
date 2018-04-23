\name{Isomap}
\alias{Isomap}
\title{
	Isomap
}
\description{
	Computes the Isomap embedding as introduced in 2000 by Tenenbaum, de Silva and Langford.
}
\usage{
Isomap(data, dims = 2, k, mod = FALSE, plotResiduals = FALSE, verbose = TRUE)
}
\arguments{
  \item{data}{
	N x D matrix (N samples, D features)
}
  \item{dims}{
	vector containing the target space dimension(s)
}
  \item{k}{
	number of neighbours
}
  \item{mod}{
	use modified Isomap algorithm
}
  \item{plotResiduals}{
	show a plot with the residuals between the high and the low dimensional data
}
  \item{verbose}{
	show a summary of the embedding procedure at the end
}
}
\details{
  Isomap is a nonlinear dimension reduction technique, that preserves
  global properties of the data. That means, that geodesic distances
  between all samples are captured best in the low dimensional
  embedding. \cr
  This R version is based on the Matlab implementation by Tenenbaum and
  uses Floyd's Algorithm to compute the neighbourhood graph of shortest
  distances, when calculating the geodesic distances. \cr
  A modified version of the original Isomap algorithm is included. It
  respects nearest and farthest neighbours. \cr
  To estimate the intrinsic dimension of the data, the function can plot
  the residuals between the high and the low dimensional data for a
  given range of dimensions.
}
\value{
	It returns a N x dim matrix (N samples, dim features) with the reduced input data (list of several matrices if more than one dimension was specified)
}
\references{
	Tenenbaum, J. B. and de Silva, V. and Langford, J. C., "A global geometric framework for nonlinear dimensionality reduction.", 2000;
	Matlab code is available at http://waldron.stanford.edu/~isomap/
}
\author{
	Christoph Bartenhagen
}