\name{stochasticProfilingML}
\alias{stochasticProfilingML}
\title{
User prompt for maximum likelihood estimation of stochastic profiling model
}
\description{
Maximum likelihood estimation for the parameters in the stochastic profiling model. The user is prompted to input the data and all required settings. He or she hence does not have to delve into the structure of this package.
}
\usage{
stochasticProfilingML()
}
\details{
The function performs maximum likelihood estimation for the parameters in the stochastic profiling model. The user is prompted to input the data and all required settings, i.e. especially the data, but also the exact model (LN-LN, rLN-LN or EXP-LN), the number of populations, the number of cells per sample etc. The data can either be entered manually, or it can be read from a file, or the user enters the name of a variable which contains the data.
}
\value{
A list as returned by \code{stochprof.loop}, i.e. the following components:
\item{mle}{maximum likelihood estimate}
\item{neg-loglikeli}{value of the negative log-likelihood function at maximum likelihood estimate}
\item{ci}{approximate marginal maximum likelihood confidence intervals for the maximum likelihood estimate}
\item{pargrid}{matrix containing parameter combinations and according values of the target function}
\item{bic}{Bayesian information criterion value}
\item{adj.bic}{adjusted Bayesian information criterion value which takes into account the numbers of parameters that were estimated during the preanalysis of a gene cluster. Is only calculated if parameter \code{subgroups} is given, otherwise set to NULL.}
\item{pen}{penalization for densities not fulfilling required constraints. If \code{use.constraints} is FALSE, this has no practical meaning. If \code{use.constraints} is TRUE, this value is included in \code{loglikeli}.}
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ stochastic profiling }
\keyword{ maximum likelihood estimation }
\keyword{ user prompt }
