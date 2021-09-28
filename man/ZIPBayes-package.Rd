\name{ZIPBayes-package}
\alias{ZIPBayes-package}
\alias{ZIPBayes}
\docType{package}
\title{
Bayesian Methods in Zero-inflated Poisson Model
}
\description{
  Implementation of zero-inflated Poisson (ZIP) model in Bayesian methods with data augmentation strategy. The R package is general to different scenarios, including an ordinary scenario for zero-inflated Poisson Model, the scenario with measurement error in response, the scenario with internal or external validation data available.
}
\details{
This package implemented the zero-inflated Poisson model based on a Bayesian framework. The method is implemented with a Monte Carlo Markov Chain (MCMC) approach with a data augmentation strategy. The package is integrated with C++ to improve the computing speed. It mainly contains four main functions. The function ZIPBayes corresponds to ordinary zero-inflated count data and no measurement error is considered. The function ZIPBayes_MErr considers the case where the response is subject to measurement error as a model by Qihuang Zhang (2020). The function ZIPBayes_Int and ZIPBayes_Ext are corresponding to the case where the internal or external validation data are available, respectively. Other helper functions are also contained in this packages, such as summarizing the trace from the MCMC algorithm, plotting the trace plot, etc.
}
\author{
Qihuang Zhang and Grace Y. Yi.

Maintainer: Qihuang Zhang <qihuang.zhang@uwaterloo.ca>
}

\keyword{package}
\keyword{Bayesian Methods}
\keyword{zero-inflated Poisson model}
\keyword{measurement error model}

