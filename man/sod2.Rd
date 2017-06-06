\name{sod2}
\alias{sod2}
\docType{data}
\title{
Measurements from the detoxifying enzyme, SOD2
}
\description{
Real 10-cell samplings from the detoxifying enzyme, SOD2. The dataset contains the measurements of SOD2 expression by qPCR in 81 random samplings of 10 ECM-attached cells.
}
\usage{data(sod2)}
\format{
  The format is:
 num [1:81] 0.603 0.873 0.204 1 3.001 ...
}
\source{
"Parameterizing cell-to-cell regulatory heterogeneities via stochastic transcriptional profiles" by Sameer S Bajikar*, Christiane Fuchs*, Andreas Roller, Fabian J Theis^ and Kevin A Janes^: PNAS 2014, 111(5), E626-635 (* joint first authors, ^ joint last authors)
}
\references{
"Parameterizing cell-to-cell regulatory heterogeneities via stochastic transcriptional profiles" by Sameer S Bajikar*, Christiane Fuchs*, Andreas Roller, Fabian J Theis^ and Kevin A Janes^: PNAS 2014, 111(5), E626-635 (* joint first authors, ^ joint last authors)
}
\examples{
data(sod2)
hist(sod2,breaks=seq(0,7,0.5),col="grey")
}
\keyword{ datasets }
\keyword{ SOD2 }
\keyword{ stochastic profiling }