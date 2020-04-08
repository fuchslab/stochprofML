\name{mix.d.sum.of.mixtures.rLNLN}
\alias{mix.d.sum.of.mixtures.rLNLN}
\title{
Density of the sum of mixtures of lognormal random variables weighted by all possible summands
}
\description{
Density and random generation of a sum of i.i.d. random variables, where each random variable is from the following mixture distribution: With probability p_i, it is of type i. In that case, it is lognormally distributed with log-mean mu_i and log-standard deviation sigma_i.  The density is somehow a "mixed" one, as for all values in \code{n} the density of the random variable is calculated and the weighted average is taken to be the density of this specific value.
}
\usage{
mix.d.sum.of.mixtures.rLNLN(y, n.vector, p.vector, mu.vector, sigma.vector)
}
\arguments{
  \item{y}{
the argument at which the density is evaluated
}
  \item{n.vector}{
the number of random variables entering each sum (in the considered application: number of cells per tissue sample). This can also be a vector stating how many cells are in each sample separatly
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
}
\details{
The lengths of p.vector, mu.vector and sigma.vector have to be identical. Their lengths automatically determine the number of different types.
}
\value{
'mix.d.sum.of.mixtures.LNLN' gives the density of a random variable originating from one of the tissue samples in the mixed n-vector.
}
\references{
"Parameterizing cell-to-cell regulatory heterogeneities via stochastic transcriptional profiles" by Sameer S Bajikar*, Christiane Fuchs*, Andreas Roller, Fabian J Theis^ and Kevin A Janes^: PNAS 2014, 111(5), E626-635 (* joint first authors, ^ joint last authors) <doi:10.1073/pnas.1311647111>

"Pheno-seq - linking visual features and gene expression in 3D cell culture systems" by Stephan M. Tirier, Jeongbin Park, Friedrich Preusser, Lisa Amrhein, Zuguang Gu, Simon Steiger, Jan-Philipp Mallm, Teresa Krieger, Marcel Waschow, Bjoern Eismann, Marta Gut, Ivo G. Gut, Karsten Rippe, Matthias Schlesner, Fabian Theis, Christiane Fuchs, Claudia R. Ball, Hanno Glimm, Roland Eils & Christian Conrad: Sci Rep 9, 12367 (2019) <doi:10.1038/s41598-019-48771-4>
}
\author{
Lisa Amrhein, Christiane Fuchs

Maintainer: Lisa Amrhein <amrheinlisa@gmail.com>
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ probability density function }
\keyword{ mixture of lognormals }
\keyword{ sum of lognormals }
\keyword{ stochastic profiling }
