\name{datasimInt}
\docType{data}
\alias{datasimInt}
\title{Toy example data - main study and internal validation study}
\description{
  This data set gives an example data for the illustration of usage of \code{\link{ZIPInt}} function. The dataset contains a list of main data and internal validation data. 
}
\usage{data(datasimInt)
}
\format{A list of two data.frames. The first data.frame, named ``main'', corresponds to the main data with 6 columns. Same as the \code{\link{datasim}}, ``Ystar'' refers to the error-prone response. ``Y'' refers to the true count response. ``X1'' and ``X2'' are covariates in the response model. ``Zplus'' and ``Zplus'' are the covariates for the measurement error model. The second data.frame corresponds to the validation data with 7 columns. }
\keyword{datasets}