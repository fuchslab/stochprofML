\name{d.sum.of.mixtures.EXPLN}
\alias{d.sum.of.mixtures.EXPLN}
\alias{r.sum.of.mixtures.EXPLN}
\title{
Sums of mixtures of zero, one or more lognormal random variables and one exponential random variable
}
\description{
Density and random generation of a sum of i.i.d. random variables, where each random variable is from the following mixture distribution: With probability p_i, it is of type i. For all but the largest i, it is lognormally distributed with log-mean mu_i and log-standard deviation sigma_i. Otherwise it is exponentially distributed with rate lambda.
}
\usage{
d.sum.of.mixtures.EXPLN(y, n, p.vector, mu.vector, sigma.vector, lambda,
   logdens = T)
r.sum.of.mixtures.EXPLN(k, n, p.vector, mu.vector, sigma.vector, lambda, N.matrix)
}
\arguments{
  \item{y}{
the argument at which the density is evaluated
}
  \item{k}{
number of i.i.d. random variables returned by this function (in the considered application: number of tissue samples)
}
  \item{n}{
the number of random variables entering each sum (in the considered application: number of cells per tissue sample). This can also be a vector stating how many cells are in each sample seperatly.
}
  \item{p.vector}{
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
  \item{N.matrix}{
optional. Matrix, that shows the decomposition of samples to be generated. Each row stands for a future sample and the columns show the decomposition of types for the specific sample such that the row sum should be n (either a value, i.e. all the same or a vector).
}
}
\details{
The lengths of mu.vector and sigma.vector have to be identical. p.vector has to have one component more. Its length automatically determines the number of different types. lambda has to be a scalar.
}
\value{
'd.sum.of.mixtures.EXPLN' gives the density, and 'r.sum.of.mixtures.EXPLN' generates random variables.
}
\references{
"Parameterizing cell-to-cell regulatory heterogeneities via stochastic transcriptional profiles" by Sameer S Bajikar*, Christiane Fuchs*, Andreas Roller, Fabian J Theis^ and Kevin A Janes^: PNAS 2014, 111(5), E626-635 (* joint first authors, ^ joint last authors)
}
\author{
Lisa Amrhein, Christiane Fuchs

Maintainer: Christiane Fuchs <christiane.fuchs@helmholtz-muenchen.de>
}
\examples{
# generate random variables
p <- c(0.25,0.75)
mu <- 2
sigma <- 0.3
lambda <- 5

stochprofML:::set.model.functions("EXP-LN")

r <- r.sum.of.mixtures.EXPLN(10^4,10,p,mu,sigma,lambda)
hist(r,xlab="Sum of mixtures of lognormals",freq=FALSE,breaks=100,ylim=c(0,0.2))

# plot according theoretical density function
x <- seq(round(min(r)),round(max(r)),(round(max(r))-round(min(r)))/500)
y <- d.sum.of.mixtures.EXPLN(x,10,p,mu,sigma,lambda,logdens=FALSE)
lines(x,y,col="blue",lwd=3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ random number generator }
\keyword{ probability density function }
\keyword{ mixture of lognormals }
\keyword{ sum of lognormals }
\keyword{ stochastic profiling }
