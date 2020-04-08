\name{stochprofML-package}
\alias{stochprofML-package}
\alias{stochprofML}
\docType{package}
\title{
Stochastic Profiling using Maximum Likelihood Estimation
}
\description{
This package accompanies the paper "Parameterizing cell-to-cell regulatory heterogeneities via stochastic transcriptional profiles" by Sameer S Bajikar, Christiane Fuchs, Andreas Roller, Fabian J Theis and Kevin A Janes (PNAS 2014, 111(5), E626-635). In that work, we measure expression profiles from small heterogeneous populations of cells. Each cell is assumed to be from a mixture of lognormal and exponential distributions (see details). We perform maximum likelihood estimation in order to infer the mixture ratio and the parameters of the lognormal/exponential distributions from the cumulated expression measurements.
}
\details{
\tabular{ll}{
Package: \tab stochprofML\cr
Type: \tab Package\cr
Version: \tab 2.0\cr
Date: \tab 2018-05-07\cr
License: \tab GPL (>= 2)\cr
}
There are three stochastic profiling models: The lognormal-lognormal (LN-LN) model assumes that each cell is from a mixture of one or more lognormal distributions with different log-means but identical log-standard deviations. In the relaxed lognormal-lognormal (rLN-LN) model, the log-standard deviations are not necessarily identical. The exponential-lognormal (EXP-LN) model considers the mixture of zero, one or more lognormal distributions and one exponential distribution.

The user can immediately start with data generation and model estimation using the two functions \code{stochasticProfilingData} and \code{stochasticProfilingML}, respectively. These functions prompt the user to input the settings/data. Typical analyses are then performed without the user having to delve into the structure of this package.

When not using one the above two functions, parameters can be estimated calling \code{stochprof.loop}, which again utilizes three other functions: \code{stochprof.search.LNLN}/ \code{stochprof.search.rLNLN}/ \code{stochprof.search.EXPLN} in order to calculate and locally optimize the likelihood function; \cr \code{stochprof.results.LNLN}/ \code{stochprof.results.rLNLN}/ \code{stochprof.results.EXPLN} for evaluating these results; and \code{calculate.ci.LNLN}/ \code{calculate.ci.rLNLN}/ \code{calculate.ci.EXPLN} for calculating confidence intervals.

Two essential functions are \code{r.sum.of.mixtures.LNLN}/ \code{r.sum.of.mixtures.rLNLN}/ \cr \code{r.sum.of.mixtures.EXPLN} and \code{d.sum.of.mixtures.LNLN}/ \code{d.sum.of.mixtures.rLNLN}/ \cr \code{d.sum.of.mixtures.EXPLN} for the density and random number generation of the distribution assumed for all measurements in the stochastic profiling model.

Version 2.0 shows one significant extension to previous versions: The number of cells \code{n} in the input data does no longer need to be the same over all samples.

The package provides four datasets: \code{sod2}, containing real measurements for one gene, and \cr \code{toycluster.LNLN}/ \code{toycluster.rLNLN}/ \code{toycluster.EXPLN}, containing artificial data for 12 genes generated with the three stochastic profiling models.

Examples for typical analyses are given below.
}
\references{
"Parameterizing cell-to-cell regulatory heterogeneities via stochastic transcriptional profiles" by Sameer S Bajikar*, Christiane Fuchs*, Andreas Roller, Fabian J Theis^ and Kevin A Janes^: PNAS 2014, 111(5), E626-635 (* joint first authors, ^ joint last authors) <doi:10.1073/pnas.1311647111>

"Pheno-seq - linking visual features and gene expression in 3D cell culture systems" by Stephan M. Tirier, Jeongbin Park, Friedrich Preusser, Lisa Amrhein, Zuguang Gu, Simon Steiger, Jan-Philipp Mallm, Teresa Krieger, Marcel Waschow, Bjoern Eismann, Marta Gut, Ivo G. Gut, Karsten Rippe, Matthias Schlesner, Fabian Theis, Christiane Fuchs, Claudia R. Ball, Hanno Glimm, Roland Eils & Christian Conrad: Sci Rep 9, 12367 (2019) <doi:10.1038/s41598-019-48771-4>
}
\author{
Lisa Amrhein, Christiane Fuchs

Maintainer: Lisa Amrhein <amrheinlisa@gmail.com>
}
\keyword{ package }
\keyword{ stochastic profiling }
\examples{
\dontrun{
# Generate synthetic data. The user is prompted to input all settings.
stochasticProfilingData()

# Estimate a stochastic profiling model. The user is prompted to input the data
# and other parameters.
stochasticProfilingML()

# Generate a synthetic dataset (without measurement error) for one gene
# and estimate the parameters from this data.
generate.toydata()

# Estimate the model parameters for the SOD2 dataset.
analyze.sod2()

# Estimate the model parameters for the 12-gene toycluster.
analyze.toycluster()
}
}
