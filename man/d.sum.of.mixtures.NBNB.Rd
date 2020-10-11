\name{d.sum.of.mixtures.NBNB}
\alias{d.sum.of.mixtures.NBNB}
\alias{r.sum.of.mixtures.NBNB}
\title{
Sums of mixtures of negative binomial random variables
}
\description{
Density and random generation of a sum of i.i.d. random variables, where each random variable is from the following mixture distribution: With probability p_i, it is of type i. In that case, it is negative binomially distributed with size alpha_i and mean mu_i.
}
\usage{
d.sum.of.mixtures.NBNB(y, n, p.vector, size.vector, mu.vector, logdens = TRUE)
r.sum.of.mixtures.NBNB(k, n, p.vector, size.vector, mu.vector, N.matrix)
}
\arguments{
  \item{y}{
the argument at which the density is evaluated
}
  \item{k}{
number of i.i.d. random variables returned by this function (in the considered application: number of tissue samples)
}
  \item{n}{
the number of random variables entering each sum (in the considered application: number of cells per tissue sample). This can also be a vector stating how many cells are in each sample separatly.
}
  \item{p.vector}{
vector (p1,p2,..,pT) containing the probabilities for each type of cell. Its elements have to sum up to one
}
  \item{size.vector}{
vector (size1,...,sizeT) containing the sizes for each type
}
  \item{mu.vector}{
vector (mu1,mu2,...,muT) containing the means for each type
}
  \item{logdens}{
if TRUE, the log of the density is returned
}
  \item{N.matrix}{
optional. Matrix, that shows the decomposition of samples to be generated. Each row stands for a future sample and the columns show the decomposition of types for the specific sample such that the row sum should be n (either a value, i.e. all the same or a vector).
}
}
\details{
The lengths of p.vector, mu.vector and size.vector have to be identical. Their lengths automatically determine the number of different types.
}
\value{
'd.sum.of.mixtures.NBNB' gives the density, and 'r.sum.of.mixtures.NBNB' generates random variables.
}
\author{
Lisa Amrhein, Christiane Fuchs

Maintainer: Lisa Amrhein <amrheinlisa@gmail.com>
}
\examples{
# generate random variables
p <- c(0.25, 0.75)
mu <- c(2, 5)
size <- c(0.3, 0.1)

set.model.functions("NB-NB")

r <- r.sum.of.mixtures.NBNB(10^4, 10, p, size.vector = size, mu.vector = mu)
hist(r, xlab = "Sum of mixtures of lognormals", freq = FALSE, breaks = 100, ylim = c(0, 0.04))

# plot according theoretical density function
x <- seq(round(min(r)), round(max(r)),(round(max(r))-round(min(r)))/500)
y <- d.sum.of.mixtures.NBNB(x, 10, p, mu, size, logdens = FALSE)
lines(x, y, col="blue", lwd = 3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ random number generator }
\keyword{ probability mass function }
\keyword{ mixture of negative binomials }
\keyword{ sum of negative binomials}
\keyword{ stochastic profiling }
