\name{ZIPExt}
\alias{ZIPExt}
\title{
Zero-inflated Poisson model under measurement error and external validation data are available
}
\description{
The function implements the MCMC algorithm with data augmentation to estimate the parameters in the zero-inflated Poisson model while correcting for the measurement error arising from the responses. The function returns the trace of the sampled parameters in each interaction. To obtain the summary estimation, use \code{summary}().
}
\usage{
ZIPExt (Ystar, Covarmainphi, Covarmainmu, Covarplus, Covarminus,
           Ystarval, Yval, Covarvalplus, Covarvalminus,
           betaphi, betamu, alphaplus, alphaminus,
           Uibound = c(7,11),
           priorgamma, priormu, priorSigma, 
           propsigmaphi,  propsigmamu = propsigmaphi, 
           propsigmaplus = propsigmaphi, propsigmaminus = propsigmaphi,
           seed = 1, nmcmc = 500)
}
\arguments{
  \item{Ystar}{
a count vector of length \eqn{n} specifying the error-prone response in the zero-inflated Poisson model.
}
  \item{Covarmainphi}{
a \eqn{n x p_1} dimensional data.frame or matrix of the covariate data with respect to the probability component of the zero-inflated Poisson model
}
  \item{Covarmainmu}{
a \eqn{n x p_2} dimensional data.frame or matrix of the covariate data with respect to the mean component of the zero-inflated Poisson model
}
  \item{Covarplus}{
a \eqn{n x q_1} dimensional data.frame or matrix of the covariate data for the measurement error model of the add-in error process
}
  \item{Covarminus}{
a \eqn{n x q_2} dimensional data.frame or matrix of the covariate data for the measurement error model of the leave-out error process
}
  \item{Ystarval}{
a count vector of length \eqn{m} specifying the  error-prone response in the validation data.
}
  \item{Yval}{
a count vector of length \eqn{m} specifying the precisely measured response in the validation data.
}
\item{Covarvalplus}{
a \eqn{m x q_1} dimensional data.frame or matrix of the covariate for validation data of the add-in error process
}
  \item{Covarvalminus}{
a \eqn{m x q_2} dimensional data.frame or matrix of the covariate  for validation data of the leave-out error process
}
  \item{betaphi}{
a vector of length \eqn{p_1} specifying the initial values of the parameters in the probability component of the zero-inflated Poisson model
}
  \item{betamu}{
a vector of length \eqn{p_2} specifying the initial values of the parameters in the probability component of the zero-inflated Poisson model
}
 \item{alphaplus}{
a vector of length \eqn{q_1} specifying the initial values of the parameters for the measurement error model of the add-in error process
}
  \item{alphaminus}{
a vector of length \eqn{q_2} specifying the initial values of the parameters in the probability component of the leave-out error process
}
  \item{Uibound}{
a vector of length \eqn{2} specifying the maximum number of the count in the inverse sampling method
}
  \item{priorgamma}{
a vector of length \eqn{2} specifying the two parameters of gamma prior
}
  \item{priormu}{
a vector of length \eqn{q_2} specifying the mean vector of the normal prior for the measurement error model of the leave-out error process
}
  \item{priorSigma}{
a vector of length \eqn{q_2} specifying the standard erros of the normal prior for the measurement error model of the leave-out error process
}
  \item{propsigmaphi}{
a vector of length \eqn{p_1} specifying the standard error of the Gaussian proposal distribution for the parameters corresponds to the probability component. 
}
  \item{propsigmamu}{
a vector of length \eqn{p_2} specifying the standard error of the Gaussian proposal distribution for the parameters corresponds to the mean component. 
}
 \item{propsigmaplus}{
a vector of length \eqn{q_1} specifying the standard error of the Gaussian proposal distribution for the parameters corresponds to the add-in error process. 
}
 \item{propsigmaminus}{
a vector of length \eqn{q_1} specifying the standard error of the Gaussian proposal distribution for the parameters corresponds to the leave-out error process. 
}
  \item{seed}{
a numeric value specifying the seed for random generator
}
  \item{nmcmc}{
a integer specify the number of the generation of MCMC algorithm
}
}
\details{
Comparing to the ZIPMErr function, this function has an addition component -- validation data. Here, the argument ``Ystarval", ``Yval", ``Covarvalplus", ``Covarvalminus", are new for the sceanrio with external validation.  
}
\value{
\item{BayesResults}{the list of trace of generated parameters for each component of the models. Data frame "betaphi_trace" corresponds to the probability component of ZIP response model; "betamu_trace" refers to the mean component of the ZIP response model. Data frames "alphaplus_trace" and "alphaminus_trace", respectively, correspond to the add-in error and leave-out error process in the measruement error model.}
}

\references{Zhang, Qihuang. "Inference Methods for Noisy Correlated Responses with Measurement Error." (2020).}

\author{
Qihuang Zhang and Grace Y. Yi
}


\seealso{
\code{\link[stats]{glm}}
}
\examples{
## load data
data(datasimExt)

example_ZIP_Ext <-  ZIPExt (Ystar = datasimExt$main$Ystar,
                            Covarmainphi = datasimExt$main[,c("intercept","X1")],
                            Covarmainmu = datasimExt$main[,c("intercept","X2")],
                            Covarplus = datasimExt$main[,c("intercept","Zplus")],
                            Covarminus = datasimExt$main[,c("intercept","Zminus")],
                            Ystarval = datasimExt$validation$Ystar, 
                            Yval = datasimExt$validation$Y,
                            Covarvalplus = datasimExt$validation[,3:4],
                            Covarvalminus = datasimExt$validation[,3:4],
                            betaphi = c(0.7,-0.7), betamu = c(1,-1.5), 
                            alphaplus = c(0,0), alphaminus=c(0,0),
                            priorgamma = c(0.001,0.001), priormu = c(0,0),
                            priorSigma = c(1,1), propsigmaphi = c(0.05,0.05), 
                            nmcmc = 10) 

summary(example_ZIP_Ext)
}
\keyword{Core}
\keyword{Main function}
