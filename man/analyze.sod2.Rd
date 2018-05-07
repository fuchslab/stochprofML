\name{analyze.sod2}
\alias{analyze.sod2}
\title{
Analysis of SOD2 data in stochastic profiling model
}
\description{
Estimation of the model parameters for the SOD2 dataset provided in this package.
}
\usage{
analyze.sod2(model = "LN-LN", TY = 2, use.constraints = F)
}
\arguments{
  \item{model}{
model for which one wishes to estimate the parameters: "LN-LN", "rLN-LN" or "EXP-LN"
}
  \item{TY}{
number of types of cells that is assumed in the stochastic model
}
  \item{use.constraints}{
if TRUE, constraints on the individual population densities are applied; see \code{penalty.constraint.LNLN}, \code{penalty.constraint.rLNLN} and \cr \code{penalty.constraint.EXPLN} for details.
}
}
\details{
The \code{sod2} dataset contains real 10-cell samplings from the detoxifying enzyme, SOD2. This function estimates the parameters of the stochastic profiling models for this data. At the end, it graphically represents a histogram of the SOD2 data together with the estimated probability density function.
}
\value{
A list as returned by \code{stochprof.loop}, i.e. the following components:
\item{mle}{maximum likelihood estimate}
\item{neg-loglikeli}{value of the negative log-likelihood function at maximum likelihood estimate}
\item{ci}{approximate marginal maximum likelihood confidence intervals for the maximum likelihood estimate}
\item{pargrid}{matrix containing parameter combinations and according values of the target function}
\item{bic}{Bayesian information criterion value}
\item{adj.bic}{adjusted Bayesian information criterion value which takes into account the numbers of parameters that were estimated during the preanalysis of a gene cluster (not applicable here, hence NULL)}
\item{pen}{penalization for densities not fulfilling required constraints. If \code{use.constraints} is FALSE, this has no practical meaning. If \code{use.constraints} is TRUE, this value is included in \code{loglikeli}.}
}
\references{
"Parameterizing cell-to-cell regulatory heterogeneities via stochastic transcriptional profiles" by Sameer S Bajikar*, Christiane Fuchs*, Andreas Roller, Fabian J Theis^ and Kevin A Janes^: PNAS 2014, 111(5), E626-635 (* joint first authors, ^ joint last authors)
}
\author{
Lisa Amrhein, Christiane Fuchs

Maintainer: Christiane Fuchs <christiane.fuchs@helmholtz-muenchen.de>
}
\note{
Executing this function as it is (i.e. without any parallelization) takes approx. 5-6 minutes (LN-LN model), 10 minutes (rLN-LN model) or an hour (EXP-LN model) on a standard computer.
}
\keyword{ stochastic profiling }
\keyword{ SOD2 }
