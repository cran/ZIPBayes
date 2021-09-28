\name{summary.ZIPBayes}
\alias{summary.ZIPBayes}

\title{
Summarizing the trace output from the MCMC algorithm
}
\description{
This function is a method for \code{ZIPBayes} object. It summarize the trace output from the main functions into interesting summary statistics, such as mean, median, confidence interval, and highest density region (HDR),
}
\usage{
## S3 method for class 'ZIPBayes'

\method{summary}{ZIPBayes}(object, burnin = 1, thinperiod = 1, confidence.level = 0.95, ...) 
}

\arguments{
  \item{object}{
  the ``ZIPBayes'' object gotten from the main function.}
  \item{burnin}{
  the number of records to be discarded as the early period of MCMC algorithm. Default is 1, meaning the first data point will be discarded when calculating the summary statistics.  
}
\item{thinperiod}{
  the number of period in the thining periord. The results will be picked every this number. Default is 1, meaning no thining will be done. See details.
}
\item{confidence.level}{
  the confidence level for the constructed confidence interval. Default is 0.95.
}
  \item{\dots}{
  other arguments passed to the function.
}
}
\details{
 This function summarizes the tracing results produced by \code{\link[ZIPBayes]{ZIP}}, \code{\link[ZIPBayes]{ZIPMErr}}, \code{\link[ZIPBayes]{ZIPExt}}, and \code{\link[ZIPBayes]{ZIPInt}}.
To diminish the influence of the starting values, we generally discard the first portion of each sequence and focus attention on the remaining. The argument \code{burnin} is set to control the number of steps to be discarded.

Another issue that sometimes arises, once approximate convergence has been reached, is whether to thin the sequences by keeping every \eqn{k}-th simulation draw from each sequence  and discarding the rest. The argument \code{thinperiod} is used to set \eqn{k} here.

}
\value{
\item{ZIPBayes}{a list of summary for each data set. "HDR_LB" and "HDR_UB" respectively respresents the lower and upper bound of the high density region.}
}

\author{
 Qihuang Zhang and Grace Y. Yi.
}

\seealso{
\code{\link[ZIPBayes]{ZIP}}, \code{\link[ZIPBayes]{ZIPMErr}}, \code{\link[ZIPBayes]{ZIPExt}}, \code{\link[ZIPBayes]{ZIPInt}}
}

\examples{
## Please see the example in ZIP() function
}
\keyword{Helper function}

