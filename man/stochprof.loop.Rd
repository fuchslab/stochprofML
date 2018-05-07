\name{stochprof.loop}
\alias{stochprof.loop}
\title{
Maximum likelihood estimation for the parameters in the stochastic profiling model
}
\description{
Maximum likelihood estimation for the parameters in the stochastic profiling model. Because the log-likelihood function is potentially multimodal, no straightforward use of gradient-based approaches for finding globally optimal parameter combinations is possible. To tackle this challenge, this function performs a two-step estimation procedure.
}
\usage{
stochprof.loop(model, dataset, n, TY, genenames = NULL, fix.mu = F,
   fixed.mu, par.range = NULL, prev.result = NULL, loops = 5,
   until.convergence = T, print.output = T, show.plots = T,
   plot.title = "", pdf.file, use.constraints = F, subgroups)
}
\arguments{
  \item{model}{
model for which one wishes to estimate the parameters: "LN-LN", "rLN-LN" or "EXP-LN"
}
  \item{dataset}{
matrix which contains the cumulated expression data over all cells in a tissue sample. Columns represent different genes, rows represent different tissue samples.
}
  \item{n}{
number of cells taken from each tissue sample. This can also be a vector stating how many cells are in each sample seperatly.
}
  \item{TY}{
number of types of cells that is assumed in the stochastic model
}
  \item{genenames}{
names of the genes in the dataset. For genenames==NULL, the genes will simply be enumerated according to the column numbers in the dataset.
}
  \item{fix.mu}{
if TRUE, the log-means of the lognormal distributions are kept fixed in the estimation procedure. Otherwise, they are to be estimated.
}
  \item{fixed.mu}{
vector containing the values to which the log-means should be fixed if fix.mu==T. The order of components is as follows:

(mu_type_1_gene_1, mu_type_1_gene_2, ...,

 mu_type_2_gene_1, mu_type_2_gene_2, ...).

This argument needs to be specified only when fix.mu==T.
}
  \item{par.range}{
range from which the parameter values should be randomly drawn if there is no knowledge from previous iterations of the search algorithm available. This is a matrix with the number of rows being equal to the number of model parameters. The first column contains the lower bound, the second column the upper bound. If par.range==NULL, some rather large range is defined.
}
  \item{prev.result}{
can contain results from former calls of this function
}
  \item{loops}{
maximal number of loops carried out in the estimation procedure. Each loops involves various methods to determine the high-likelihood region.
}
  \item{until.convergence}{
if TRUE, the estimation process is terminated if there had been no improvement concerning the value of the target function between two consecutive loops. Otherwise, the estimation procedure is terminated according to the parameter "loops".
}
  \item{print.output}{
if TRUE, interim results of the grid search and numerical optimization are printed into the console throughout the estimation procedure
}
  \item{show.plots}{
if TRUE, profile log-likelihood plots are displayed at the end of the estimation procedure
}
  \item{plot.title}{
title of each plot if show.plots==T
}
  \item{pdf.file}{
optional filename. If this is not missing and show.plots==T, the profile log-likelihoods will be plotted into this file.
}
  \item{use.constraints}{
if TRUE, constraints on the individual population densities are applied; see \code{penalty.constraint.LNLN}, \code{penalty.constraint.rLNLN} and \cr \code{penalty.constraint.EXPLN} for details.
}
  \item{subgroups}{
list of sets of gene numbers. This parameter should be given only when the present call of \code{stochprof.loop} is based on a subanalysis of the subgroups of genes with non-fixed mu. The parameter is used only for calculation of the adjusted BIC which takes into account the number of parameters that had to be estimated during the whole estimation procedure: First, for each of the subclusters, and then for the final analysis.
}
}
\details{
This function carries out maximum likelihood estimation for the parameters of the stochastic profiling model. Because the log-likelihood function is potentially multimodal, no straightforward use of gradient-based approaches for finding globally optimal parameter combinations is possible. To tackle this challenge, this function performs a two-step estimation procedure: First, it computes the log-likelihood function at randomly drawn parameter combinations to identify high-likelihood regions in parameter space at computationally  low cost. Then, it uses the Nelder-Mead algorithm to identify local maxima of the likelihood function. The starting values for this algorithm are randomly drawn from the high-likelihood regions determined in the first step. To further localize the global optimum, the function again performs grid searches of the parameter space, this time around the optimum identified by the Nelder-Mead algorithm. This search creates another space to identify high-likelihood regions, which are then used to seed another Nelder-Mead optimization.
}
\value{
A list with the following components:
\item{mle}{maximum likelihood estimate}
\item{neg-loglikeli}{value of the negative log-likelihood function at maximum likelihood estimate}
\item{ci}{approximate marginal maximum likelihood confidence intervals for the maximum likelihood estimate}
\item{pargrid}{matrix containing parameter combinations and according values of the target function}
\item{bic}{Bayesian information criterion value}
\item{adj.bic}{adjusted Bayesian information criterion value which takes into account the numbers of parameters that were estimated during the preanalysis of a gene cluster. Is only calculated if parameter \code{subgroups} is given, otherwise set to NULL.}
\item{pen}{penalization for densities not fulfilling required constraints. If \code{use.constraints} is FALSE, this has no practical meaning. If \code{use.constraints} is TRUE, this value is included in \code{loglikeli}.}
}
\references{
"Parameterizing cell-to-cell regulatory heterogeneities via stochastic transcriptional profiles" by Sameer S Bajikar*, Christiane Fuchs*, Andreas Roller, Fabian J Theis^ and Kevin A Janes^: PNAS 2014, 111(5), E626-635 (* joint first authors, ^ joint last authors)
}
\author{
Lisa Amrhein, Christiane Fuchs

Maintainer: Christiane Fuchs <christiane.fuchs@helmholtz-muenchen.de>
}
\keyword{ maximum likelihood estimation }
\keyword{ stochastic profiling }
