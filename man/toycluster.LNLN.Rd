\name{toycluster.LNLN}
\alias{toycluster.LNLN}
\docType{data}
\title{
Synthetic data from the LN-LN model
}
\description{
A matrix containing synthetic measurements from the stochastic profiling LN-LN model. There is data for 12 genes (columns) and 16 tissue samples (rows). Each measurement is the sum of 10 i.i.d. random variables from a mixture of lognormal distributions.
}
\usage{data(toycluster.LNLN)}
\format{
  The format is:
 num [1:16, 1:12] 0.789 4.698 4.643 8.734 12.458 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:16] "V1" "V2" "V3" "V4" ...
  ..$ : chr [1:12] "gene 1" "gene 2" "gene 3" "gene 4" ...
}
\details{
The true underlying parameters are:

TY = 2, i.e. there are two types of cells

p = (0.225, 0.775), that is the probability for cell type I and II, respectively

mu1 = (1.8853, 2.2758, 0.4748, 0.2658, 1.5745, 2.3938, 1.7389, 2.2148, 0.2104, 2.1032, 0.0638, 1.8109), that is the log-mean for cell type I for genes 1 to 12

mu2 = (-2.6637, -0.6590, -1.6308, -2.0753, -1.5786, -0.8131, -2.4872, -3.4486, -3.4865, -2.1848, -1.3868, -2.8238), that is the log-mean for cell type II for genes 1 to 12
  
sigma = 0.225, that is the log-standard deviation for both cell types
}
\source{
"Parameterizing cell-to-cell regulatory heterogeneities via stochastic transcriptional profiles" by Sameer S Bajikar*, Christiane Fuchs*, Andreas Roller, Fabian J Theis^ and Kevin A Janes^: PNAS 2014, 111(5), E626-635 (* joint first authors, ^ joint last authors)
}
\references{
"Parameterizing cell-to-cell regulatory heterogeneities via stochastic transcriptional profiles" by Sameer S Bajikar*, Christiane Fuchs*, Andreas Roller, Fabian J Theis^ and Kevin A Janes^: PNAS 2014, 111(5), E626-635 (* joint first authors, ^ joint last authors)
}
\examples{
data(toycluster.LNLN)
par(mfrow=c(3,4))
for (i in 1:ncol(toycluster.LNLN)) {
   hist(toycluster.LNLN[,i],xlab="synthetic data from LN-LN model",
      main=colnames(toycluster.LNLN)[i],col="grey")
}
par(mfrow=c(1,1))
}
\keyword{ datasets }
\keyword{ synthetic data }
\keyword{ stochastic profiling }