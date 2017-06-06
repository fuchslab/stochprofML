\name{toycluster.rLNLN}
\alias{toycluster.rLNLN}
\docType{data}
\title{
Synthetic data from the rLN-LN model
}
\description{
A matrix containing synthetic measurements from the stochastic profiling rLN-LN model. There is data for 12 genes (columns) and 16 tissue samples (rows). Each measurement is the sum of 10 i.i.d. random variables from a mixture of lognormal distributions.
}
\usage{data(toycluster.rLNLN)}
\format{
  The format is:
 num [1:16, 1:12] 3.46 2.34 3.98 3.42 3.43 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:16] "V1" "V2" "V3" "V4" ...
  ..$ : chr [1:12] "gene 1" "gene 2" "gene 3" "gene 4" ...
}
\details{
The true underlying parameters are:

TY = 2, i.e. there are two types of cells

p = (0.225, 0.775), that is the probability for cell type I and II, respectively

mu1 = (0.1287, 1.6249, 1.0075, 0.5521, 0.1200, 1.1661, 1.4261, 1.8238, 2.4261, 1.2568, 0.9342, 1.8876), that is the log-mean for cell type I for genes 1 to 12

mu2 = (-2.2181, -1.6432, -0.9966, -3.1968, -1.9852, -1.0545, -2.3596, -3.0939, -1.3195, -3.2041, -1.2185, -1.3895), that is the log-mean for cell type II for genes 1 to 12
  
sigma = (0.225,0.625), that are the log-standard deviations for the two cell types
}
\source{
"Parameterizing cell-to-cell regulatory heterogeneities via stochastic transcriptional profiles" by Sameer S Bajikar*, Christiane Fuchs*, Andreas Roller, Fabian J Theis^ and Kevin A Janes^: PNAS 2014, 111(5), E626-635 (* joint first authors, ^ joint last authors)
}
\references{
"Parameterizing cell-to-cell regulatory heterogeneities via stochastic transcriptional profiles" by Sameer S Bajikar*, Christiane Fuchs*, Andreas Roller, Fabian J Theis^ and Kevin A Janes^: PNAS 2014, 111(5), E626-635 (* joint first authors, ^ joint last authors)
}
\examples{
data(toycluster.rLNLN)
par(mfrow=c(3,4))
for (i in 1:ncol(toycluster.rLNLN)) {
   hist(toycluster.rLNLN[,i],xlab="synthetic data from rLN-LN model",
      main=colnames(toycluster.rLNLN)[i],col="grey")
}
par(mfrow=c(1,1))
}
\keyword{ datasets }
\keyword{ synthetic data }
\keyword{ stochastic profiling }