\name{penalty.constraint.NBNB}
\alias{penalty.constraint.NBNB}
\title{
Penalization for population densities that do not fulfil certain constraints for the NB-NB model
}
\description{
In order to force the individual populations to be sufficiently distinct from each other, one can perform penalized optimization. To this end, constraints on the densities are introduced (see details). If the constraints are not fulfilled, a penalization term is added to the negative log-likelihood (which is to be minimized).
}
\usage{
penalty.constraint.NBNB(dataset, parameter, smoothingpar = 10^5)
}
\arguments{
  \item{dataset}{
matrix which contains the cumulated expression data over all cells in a tissue sample. Columns represent different genes, rows represent different tissue samples.
}
  \item{parameter}{
parameter for which the penalization term is calculated. This is a vector containing p, mu and size.
}
  \item{smoothingpar}{
weight with which the penalization term is multiplied.
}
}
\details{
The constraints are as follows: There are TY densities for the TY distinct populations. For each i=1,...,(TY-1), one considers the density of population i (the higher regulatory state) and the density of population i+1 (the lower regulatory state). The density of the higher regulatory state is constrained to be greater than the density of the lower regulatory state in the domain between the mode of the high state and the largest observation in the dataset.

Introduction of this penalization term does not mean that the constraints will automatically be fulfilled. The parameter estimate will be a trade-off between a maximizer of the unconstrained likelihood function and a minimizer of the penalization function. The higher the parameter \code{smoothingpar}, the more importance is on fulfilling the constraints.
}
\value{
The population densities are compared on the above described domains. Wherever the constraint is not fulfilled, the difference between the larger and the lower density is calculated. The squares of all such differences are summed up and multiplied with \code{smoothingpar}. This value is returned.
}
\author{
Lisa Amrhein, Christiane Fuchs

Maintainer: Lisa Amrhein <amrheinlisa@gmail.com>
}
\keyword{ stochastic profiling }
\keyword{ constraints }
\keyword{ penalization }
