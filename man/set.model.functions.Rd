\name{set.model.functions}
\alias{set.model.functions}
\title{
Defines some global model-dependent functions
}
\description{
Defines which model to use in the following and sets depending on this decistion some global model-dependent functions. This function is unly needed when using lower leverl functions of the package such as d.sum.of.mixtures which does not have its own model selection input.
}
\usage{
set.model.functions(model)
}
\arguments{
  \item{model}{
model for which one wishes to estimate the parameters: "LN-LN", "rLN-LN", "EXP-LN" or "NB-NB"
}
}
\value{
Loads all model-specific functions and directs the general function to the model specific function, wich is the same but with the endings .LNLN, .rLNLN, .EXPLN. or .NBNB
}

\author{
Lisa Amrhein, Christiane Fuchs

Maintainer: Lisa Amrhein <amrheinlisa@gmail.com>
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ stochastic profiling }
\keyword{ model selection }
\keyword{ model dependent function definition }
