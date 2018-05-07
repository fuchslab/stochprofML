\name{mix.d.sum.of.mixtures.EXPLN}
\alias{mix.d.sum.of.mixtures.EXPLN}
\title{
Density of the sum of mixtures of zero, one or more lognormal random variables and one exponential random variable weighted by all possible summands
}
\description{
Density of a sum of i.i.d. random variables, where each random variable is from the following mixture distribution: With probability p_i, it is of type i. For all but the largest i, it is lognormally distributed with log-mean mu_i and log-standard deviation sigma_i. Otherwise it is exponentially distributed with rate lambda. The density is somehow a "mixed" one, as for all values in \code{n} the density of the random varibale is calculated and the weighted average is taken to be the density of this specific value.
}
\usage{
mix.d.sum.of.mixtures.EXPLN(y, n, p.vector, mu.vector, sigma.vector, lambda,
   logdens = T)
}
\arguments{
  \item{y}{
the argument at which the density is evaluated
}
  \item{k}{
number of i.i.d. random variables returned by this function (in the considered application: number of tissue samples)
}
  \item{n}{
the number of random variables entering each sum (in the considered application: number of cells per tissue sample). This can also be a vector stating how many cells are in each sample seperatly
}
  \item{p.vector}{
vector (p1,p2,..,pT) containing the probabilities for each type of cell. Its elements have to sum up to one
}
  \item{mu.vector}{
vector (mu1,mu2,...,mu(T-1)) containing the log-means for each lognormal type (types 1 to T-1)
}
  \item{sigma.vector}{
vector (sigma1,...,sigma(T-1)) containing the log-standard deviations sigma for each lognormal type (types 1 to T-1)
}
  \item{lambda}{
the rate for the exponential type (type T)
}
  \item{logdens}{
if TRUE, the log of the density is returned
}
}
\details{
The lengths of mu.vector and sigma.vector have to be identical. p.vector has to have one component more. Its length automatically determines the number of different types. lambda has to be a scalar.
}
\value{
'mix.d.sum.of.mixtures.EXPLN' gives the density of a random variable originating from one of the tissue samples in the mixed n-vector.
}

\author{
Lisa Amrhein, Christiane Fuchs

Maintainer: Christiane Fuchs <christiane.fuchs@helmholtz-muenchen.de>
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ probability density function }
\keyword{ mixture of lognormals }
\keyword{ sum of lognormals }
\keyword{ stochastic profiling }
