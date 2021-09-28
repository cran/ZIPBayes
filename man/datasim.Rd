\name{datasim}
\docType{data}
\alias{datasim}
\title{Toy example data - main study only}
\description{
  This data set gives an example data for the illustration of usage of \code{\link{ZIP}} and \code{\link{ZIPMErr}} function. The dataset contains naivedata and the design matrices for the response model, measurement error model.  
}
\usage{data(datasim)
}
\format{A data.frame of 6 columns. ``Ystar" refers to the error-prone response. ``Y" refers to the true count response. "X1" and "X2" are covariates in the response model. ``Xplus" and ``Xminus" are the covariates for the measurement error model. }
\keyword{datasets}