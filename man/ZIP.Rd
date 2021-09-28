\name{ZIP}
\alias{ZIP}
\title{
zero-inflated Poisson model
}
\description{
The function implements the MCMC algorithm with data augmentation to estimate the parameters in the zero-inflated Poisson model. The function returns the trace of the sampled parameters in each interaction. To obtain the summary estimation, use \code{summary}().
}
\usage{
ZIP(Y, Covarmainphi, Covarmainmu, 
                         betaphi, betamu, 
                         priorgamma, 
                         propsigmaphi,  propsigmamu = propsigmaphi,
                         seed = 1, nmcmc = 500)
}
\arguments{
  \item{Y}{
a count vector of length \eqn{n} specifying response in the zero-inflated Poisson model.
}
  \item{Covarmainphi}{
a \eqn{n x p_1} dimensional data.frame or matrix of  data with respect to the probability component of the zero-inflated Poisson model. 
}
  \item{Covarmainmu}{
a \eqn{n x p_2} dimensional data.frame or matrix of  data with respect to the mean component of the zero-inflated Poisson model. 
}
  \item{betaphi}{
a vector of length \eqn{p_1} specifying the initial values of the parameters in the probability component of the zero-inflated Poisson model
}
  \item{betamu}{
a vector of length \eqn{p_2} specifying the initial values of the parameters in the probability component of the zero-inflated Poisson model
}
  \item{priorgamma}{
a vector of length \eqn{2} specifying the two parameters of gamma prior
}
  \item{propsigmaphi}{
a vector of length \eqn{p_1} specifying the standard error of the Gaussian proposal distribution for the parameters corresponds to the probability component. 
}
  \item{propsigmamu}{
a vector of length \eqn{p_2} specifying the standard error of the Gaussian proposal distribution for the parameters corresponds to the mean component. 
}
  \item{seed}{
a numeric value specifying the seed for random generator
}
  \item{nmcmc}{
a integer specify the number of the generation of MCMC algorithm
}
}
\details{
The zero-inflated Poisson model involves two components, the probability components and the mean compoenents (Zhang, 2020). Argument \code{Covarmainphi}, \code{betaphi}, \code{propsigmaphi} correspond to the probability compoenent; \code{Covarmainmu}, \code{betamu}, \code{propsigmamu} correspond to the mean compoenent.
}
\value{
\item{BayesResults}{the list of trace of generated parameters for each component of the models. Data.frame "betaphi_trace" corresponds to the probability component of ZIP response model; "betamu_trace" refers to the mean component of the ZIP response model. }
}

\references{Zhang, Qihuang. "Inference Methods for Noisy Correlated Responses with Measurement Error." (2020).}

\author{
Qihuang Zhang and Grace Y. Yi
}


\seealso{
\code{\link[stats]{glm}}
}
\examples{
data(datasim)

example_ZIP <- ZIP( Y = datasim$Ystar,
                         Covarmainphi = datasim[,c("intercept","X1")],
                         Covarmainmu = datasim[,c("intercept","X2")],
                         betaphi = c(-0.7,0.7), betamu = c(1,-0.5),
                         priorgamma = rep(1,1), propsigmaphi = c(0.05,0.05),
                         nmcmc = 100)

summary(example_ZIP)
}
\keyword{Core}
\keyword{Main function}
