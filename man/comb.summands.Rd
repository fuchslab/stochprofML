\name{comb.summands}
\alias{comb.summands}
\title{
Combinations of fixed number of summands with pre-defined sum.
}
\description{
Returns all combinations of k numbers between 0 and n whose sum equals n.
}
\usage{
comb.summands(n, k)
}
\arguments{
  \item{n}{
the sum of the k summands
}
  \item{k}{
the number of summands
}
}
\details{
Returns all combinations of k numbers (non-negative integers) between 0 and n whose sum equals n. The order of the summands matters,
i.e. 2+3=5 and 3+2=5 would count as two different combinations.
}
\value{
A matrix with k columns. Each row contains a different combination of k non-negative integers which sum up to n.
}
\references{
"Parameterizing cell-to-cell regulatory heterogeneities via stochastic transcriptional profiles" by Sameer S Bajikar*, Christiane Fuchs*, Andreas Roller, Fabian J Theis^ and Kevin A Janes^: PNAS 2014, 111(5), E626-635 (* joint first authors, ^ joint last authors)
}
\author{
Christoph Kurz
}
\keyword{ combination of type numbers }
