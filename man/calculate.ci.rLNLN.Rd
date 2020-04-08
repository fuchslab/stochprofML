\name{calculate.ci.rLNLN}
\alias{calculate.ci.rLNLN}
\title{
Maximum likelihood confidence intervals for rLN-LN model
}
\description{
Calculates approximate marginal maximum likelihood confidence intervals with
significance level alpha for all parameters in the rLN-LN model.
}
\usage{
calculate.ci.rLNLN(alpha, parameter, prev.result, dataset, n, TY,
   fix.mu = F, fixed.mu)
}
\arguments{
  \item{alpha}{
the significance level
}
  \item{parameter}{
the maximum likelihood estimate around which the confidence interval is centered; if this value is missing, it is determined from prev.result. This parameter has to be on the original scale, not on the logit-/log-transformed scale as used during the estimation procedure. It has to be ((m+1)*TY-1)-dimensional, even for fix.mu==T.
}
  \item{prev.result}{
a list of previous results as returned by stochprof.results. This variable is only used when 'parameter' is missing.
}
  \item{dataset}{
a matrix which contains the cumulated expression data over all cells in a tissue sample. Columns represent different genes, rows represent different tissue samples.
}
  \item{n}{
number of cells taken from each tissue sample. This can also be a vector stating how many cells are in each sample separatly.
}
  \item{TY}{
the number of types of cells that is assumed in the stochastic model
}
  \item{fix.mu}{
if TRUE, the log-means are kept fixed in the estimation procedure. Otherwise, they are to be estimated.
}
  \item{fixed.mu}{
a vector containing the values to which the log-means should be fixed if \cr fix.mu==T. The order of components is as follows:

(mu_type_1_gene_1, mu_type_1_gene_2, ...,

 mu_type_2_gene_1, mu_type_2_gene_2, ...).

The parameter needs to be specified only when fix.mu==T.
}
}
\details{
The intervals are approximate because the function uses the formula

    [ theta_i +/- q_(1-alpha/2) * sqrt(H_ii) ],

where theta_i is the estimate of the i.th parameter, q_(1-alpha/2) is the 1-alpha/2 quantile of the standard normal distribution, H is the inverse Hessian of the negative log likelihood function evaluated at the maximum likelihood estimate; H_ii is the i.th diagonal element of H. This approximation implicitly assumes that the log likelihood function is unimodal. The confidence interval is first calculated on the transformed, unrestricted parameter space and then back-transformed to the original one.
}
\value{
Approximate marginal maximum likelihood confidence intervals for all parameter components: Each row corresponds to one parameter (in the same order as always used in the stochastic profiling model). The first column contains lower bounds, the second column upper bounds.
}
\references{
"Parameterizing cell-to-cell regulatory heterogeneities via stochastic transcriptional profiles" by Sameer S Bajikar*, Christiane Fuchs*, Andreas Roller, Fabian J Theis^ and Kevin A Janes^: PNAS 2014, 111(5), E626-635 (* joint first authors, ^ joint last authors) <doi:10.1073/pnas.1311647111>

"Pheno-seq - linking visual features and gene expression in 3D cell culture systems" by Stephan M. Tirier, Jeongbin Park, Friedrich Preusser, Lisa Amrhein, Zuguang Gu, Simon Steiger, Jan-Philipp Mallm, Teresa Krieger, Marcel Waschow, Bjoern Eismann, Marta Gut, Ivo G. Gut, Karsten Rippe, Matthias Schlesner, Fabian Theis, Christiane Fuchs, Claudia R. Ball, Hanno Glimm, Roland Eils & Christian Conrad: Sci Rep 9, 12367 (2019) <doi:10.1038/s41598-019-48771-4>
}
\author{
Lisa Amrhein, Christiane Fuchs

Maintainer: Lisa Amrhein <amrheinlisa@gmail.com>
}
\keyword{ maximum likelihood confidence interval }
\keyword{ stochastic profiling }
