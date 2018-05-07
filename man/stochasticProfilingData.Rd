\name{stochasticProfilingData}
\alias{stochasticProfilingData}
\title{
User prompt for generation and visualization of synthetic data in stochastic profiling model
}
\description{
Generation of a dataset of i.i.d measurements as considered in the stochastic profiling model. The user is prompted to input all required settings. He or she hence does not have to delve into the structure of this package.
}
\usage{
stochasticProfilingData()
}
\details{
This function generates a dataset of i.i.d synthetic measurements as considered in the stochastic profiling model. The user is prompted to input all required settings, e.g. the exact model (LN-LN, rLN-LN or EXP-LN), the number of samples, the number of populations etc. The dataset is returned as a matrix and visualized in histograms together with the data generating probability density function (pdf).
}
\value{
A matrix containing the synthetic data. Columns stand for different genes, rows for samples.
}
\references{
"Parameterizing cell-to-cell regulatory heterogeneities via stochastic transcriptional profiles" by Sameer S Bajikar*, Christiane Fuchs*, Andreas Roller, Fabian J Theis^ and Kevin A Janes^: PNAS 2014, 111(5), E626-635 (* joint first authors, ^ joint last authors)
}
\author{
Lisa Amrhein, Christiane Fuchs

Maintainer: Christiane Fuchs <christiane.fuchs@helmholtz-muenchen.de>
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ stochastic profiling }
\keyword{ data generation }
\keyword{ synthetic data }
\keyword{ user prompt }
