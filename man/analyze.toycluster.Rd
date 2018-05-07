\name{analyze.toycluster}
\alias{analyze.toycluster}
\title{
Analysis of toyclusters in stochastic profiling model
}
\description{
Estimation of the model parameters for the 12-gene toyclusters provided in this package. This is done in three steps: an optional preanalysis of the single genes, an analysis of three 4-gene subclusters, and finally the analysis of the entire 12-gene cluster.
}
\usage{
analyze.toycluster(model = "LN-LN", data.model = "LN-LN", TY = 2,
   preanalyze = T, show.plots = T, use.constraints = F)
}
\arguments{
  \item{model}{
model for which one wishes to estimate the parameters: "LN-LN", "rLN-LN" or "EXP-LN"
}
  \item{data.model}{
model which has generated the 12-gene dataset: "LN-LN", "rLN-LN" or "EXP-LN"
}
  \item{TY}{
number of types of cells that is assumed in the stochastic model
}
  \item{preanalyze}{
if TRUE, the single-gene preanalysis as described below is carried out
}
  \item{show.plots}{
if TRUE, interim results are graphically displayed. This requires the user to confirm each new plot.
}
  \item{use.constraints}{
if TRUE, constraints on the individual population densities are applied; see \code{penalty.constraint.LNLN}, \code{penalty.constraint.rLNLN} and \cr \code{penalty.constraint.EXPLN} for details.
}
}
\details{
This function carries out estimation of the model parameters for the \code{toycluster.LNLN}, \cr
\code{toycluster.rLNLN} or \code{toycluster.EXPLN} dataset. This contains perfectly observed measurements for 12 genes and 16 tissue samples, assuming 10-cell samplings and two different types of cells. The true underlying parameters are given on the help page for the datasets.

The estimation is performed in three steps:

In an optional preanalysis (carried out if preanalyze is TRUE), each gene is considered individually, i.e. for each gene the parameters are estimated (these are p, mu_1, mu_2 and sigma for LN-LN, p, mu_1, mu_2, sigma_1 and sigma_2 for rLN-LN, and p, mu, sigma and lambda for EXP-LN). This gives a rough idea about the location of the parameters at computationally low cost. This might speed up the analysis of the larger clusters. From the confidence intervals of the single-gene estimates, one can construct appropriate parameter ranges for the following step.

In the main step of the estimation procedure, the 12 genes are divided into three groups of size four. This is because the stochastic profiling model for 12 genes involves 48 (LN-LN and EXP-LN) to 49 (rLN-LN) parameters, which is computationally expensive and sometimes unreliable. Simulation studies showed that datasets comprising four genes are sufficient to estimate the log-means when there is data from 16 experiments available. For each of these 4-gene clusters, 10 (LN-LN and EXP-LN) or 11 (rLN-LN) parameters are estimated. The three groups result from a hierarchical clustering of the entire dataset. The genes numbers are (7,5,2,8), (1,3,4,10) and (9,6,12,11) for the LN-LN model, (12,9,6,11), (4,10,5,3) and (1,7,8,2) for the rLN-LN model and (11,1,10,9), (3,5,8,7) and (4,2,12,6) for the EXP-LN model.

In the final step, the log-means mu are fixed to the maximum likelihood estimates that resulted from the main step. Then there remain only p, sigma and possibly lambda to be estimated. These are inferred now.

Throughout the whole estimation process, interim results are printed into the console and, if \cr
show.plots is TRUE, graphically displayed.
}
\value{
The final result for the chosen 12-gene cluster. That is a list as returned by \code{stochprof.loop}, i.e. the following components:
\item{mle}{maximum likelihood estimate}
\item{neg-loglikeli}{value of the negative log-likelihood function at maximum likelihood estimate}
\item{ci}{approximate marginal maximum likelihood confidence intervals for the maximum likelihood estimate}
\item{pargrid}{matrix containing parameter combinations and according values of the target function}
\item{bic}{Bayesian information criterion value}
\item{adj.bic}{adjusted Bayesian information criterion value which takes into account the numbers of parameters that were estimated during the preanalysis of the gene cluster}
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
Executing this function as it is (i.e. without any parallelization) takes approx. 1-2 hours (LN-LN model), 1-3 hours (rLN-LN model) or 1 day (EXP-LN model) on a standard computer.
}
\keyword{ stochastic profiling }
\keyword{ synthetic data }
