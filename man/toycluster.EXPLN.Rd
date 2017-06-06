\name{toycluster.EXPLN}
\alias{toycluster.EXPLN}
\docType{data}
\title{
Synthetic data from the EXP-LN model
}
\description{
A matrix containing synthetic measurements from the stochastic profiling EXP-LN model. There is data for 12 genes (columns) and 16 tissue samples (rows). Each measurement is the sum of 10 i.i.d. random variables from a mixture of one lognormal and one exponential distribution.
}
\usage{data(toycluster.EXPLN)}
\format{
  The format is:
 num [1:16, 1:12] 3.77 4.87 5.05 4.45 5.35 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:16] "V1" "V2" "V3" "V4" ...
  ..$ : chr [1:12] "gene 1" "gene 2" "gene 3" "gene 4" ...
}
\details{
The true underlying parameters are:

TY = 2, i.e. there are two types of cells

p = (0.225, 0.775), that is the probability for cell type I and II, respectively

mu = (0.1223, 0.2705, 2.1457, 2.2899, 1.6791, 1.1558, 2.4035, 0.1998, 0.9648, 0.0411, 1.4798, 1.4206), that is the log-mean for cell type I for genes 1 to 12
 
sigma = 0.225, that is the log-standard deviation for type I

lambda = (5.5522, 31.5412, 21.2097, 6.1446, 49.0361, 10.9487, 29.7759, 43.8547, 35.7143, 6.5736, 24.8089, 24.7922), that is the exponential rate for cell type II for genes 1 to 12
}
\source{
"Parameterizing cell-to-cell regulatory heterogeneities via stochastic transcriptional profiles" by Sameer S Bajikar*, Christiane Fuchs*, Andreas Roller, Fabian J Theis^ and Kevin A Janes^: PNAS 2014, 111(5), E626-635 (* joint first authors, ^ joint last authors)
}
\references{
"Parameterizing cell-to-cell regulatory heterogeneities via stochastic transcriptional profiles" by Sameer S Bajikar*, Christiane Fuchs*, Andreas Roller, Fabian J Theis^ and Kevin A Janes^: PNAS 2014, 111(5), E626-635 (* joint first authors, ^ joint last authors)
}
\examples{
data(toycluster.EXPLN)
par(mfrow=c(3,4))
for (i in 1:ncol(toycluster.EXPLN)) {
   hist(toycluster.EXPLN[,i],xlab="synthetic data from EXP-LN model",
      main=colnames(toycluster.EXPLN)[i],col="grey")
}
par(mfrow=c(1,1))
}
\keyword{ datasets }
\keyword{ synthetic data }
\keyword{ stochastic profiling }