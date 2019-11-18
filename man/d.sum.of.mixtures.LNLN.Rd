\name{d.sum.of.mixtures.LNLN}
\alias{d.sum.of.mixtures.LNLN}
\alias{r.sum.of.mixtures.LNLN}
\title{
Sums of mixtures of lognormal random variables
}
\description{
Density and random generation of a sum of i.i.d. random variables, where each random variable is from the following mixture distribution: With probability p_i, it is of type i. In that case, it is lognormally distributed with log-mean mu_i and log-standard deviation sigma_i.
}
\usage{
d.sum.of.mixtures.LNLN(y, n, p.vector, mu.vector, sigma.vector, logdens = T)
r.sum.of.mixtures.LNLN(k, n, p.vector, mu.vector, sigma.vector, N.matrix)
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
  \item{mu.vector}{
vector (mu1,mu2,...,muT) containing the log-means for each type
}
  \item{sigma.vector}{
vector (sigma1,...,sigmaT) containing the log-standard deviations sigma for each type
}
  \item{logdens}{
if TRUE, the log of the density is returned
}
  \item{N.matrix}{
optional. Matrix, that shows the decomposition of samples to be generated. Each row stands for a future sample and the columns show the decomposition of types for the specific sample such that the row sum should be n (either a value, i.e. all the same or a vector).
}
}
\details{
The lengths of p.vector, mu.vector and sigma.vector have to be identical. Their lengths automatically determine the number of different types.
}
\value{
'd.sum.of.mixtures.LNLN' gives the density, and 'r.sum.of.mixtures.LNLN' generates random variables.
}
\references{
"Parameterizing cell-to-cell regulatory heterogeneities via stochastic transcriptional profiles" by Sameer S Bajikar*, Christiane Fuchs*, Andreas Roller, Fabian J Theis^ and Kevin A Janes^: PNAS 2014, 111(5), E626-635 (* joint first authors, ^ joint last authors) <doi:10.1073/pnas.1311647111>
"Pheno-seq - linking visual
  features and gene expression in 3D cell culture systems" by Stephan M. Tirier, Jeongbin Park, Friedrich Preusser, Lisa Amrhein, Zuguang Gu, Simon Steiger, Jan-Philipp Mallm, Teresa Krieger, Marcel Waschow, Bjoern Eismann, Marta Gut, Ivo G. Gut, Karsten Rippe, Matthias Schlesner, Fabian Theis, Christiane Fuchs, Claudia R. Ball, Hanno Glimm, Roland Eils & Christian Conrad: Sci Rep 9, 12367 (2019) <doi:10.1038/s41598-019-48771-4>
}
\author{
Lisa Amrhein, Christiane Fuchs

Maintainer: Lisa Amrhein <amrheinlisa@gmail.com>
}
\examples{
# generate random variables
p <- c(0.25,0.75)
mu <- c(2,-1)
sigma <- c(0.3,0.1)

stochprofML:::set.model.functions("LN-LN")

r <- r.sum.of.mixtures.LNLN(10^4,10,p,mu,sigma)
hist(r,xlab="Sum of mixtures of lognormals",freq=FALSE,breaks=100,ylim=c(0,0.2))

# plot according theoretical density function
x <- seq(round(min(r)),round(max(r)),(round(max(r))-round(min(r)))/500)
y <- d.sum.of.mixtures.LNLN(x,10,p,mu,sigma,logdens=FALSE)
lines(x,y,col="blue",lwd=3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ random number generator }
\keyword{ probability density function }
\keyword{ mixture of lognormals }
\keyword{ sum of lognormals }
\keyword{ stochastic profiling }
