\name{generate.toydata}
\alias{generate.toydata}
\title{
Generation and analysis of synthetic data in stochastic profiling model
}
\description{
Generation of a dataset of 500 i.i.d measurements as considered in the stochastic profiling model. Afterwards estimation of the model parameters and comparison of the estimates with the true value.
}
\usage{
generate.toydata(model = "LN-LN")
}
\arguments{
  \item{model}{
the chosen stochastic profiling model: "LN-LN", "rLN-LN" or "EXP-LN"
}
}
\details{
This function first generates a dataset of 500 i.i.d. 10-cell samplings as considered in the stochastic profiling models "LN-LN", "rLN-LN" and "EXP-LN". The employed parameters are TY=2 (i.e. two different types of cells are assumed) and p=c(0.2,0.8) for all models. Furthermore, mu=c(1.5,-1.5) and sigma=0.2 for the LN-LN model, mu=c(1.5,-1.5) and sigma=(0.2,0.6) for the rLN-LN model, and mu=1.5, sigma=0.2 and lambda=0.5 for the EXP-LN model. The generated data is displayed in a histogram together with the theoretical probability density function. At the end of the estimation procedure, the profile log-likelihood plots are shown. Finally, the true and the estimated probability density functions are compared and the estimation results are printed.
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
\keyword{ stochastic profiling }
\keyword{ synthetic data }
\keyword{ data generation }
