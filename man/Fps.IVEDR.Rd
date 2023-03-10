% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{Fps.IVEDR}
\alias{Fps.IVEDR}
\title{Logistic regression for observed treatment and instrument used for the (S)IVE-DR estimator.}
\usage{
Fps.IVEDR(datalist)
}
\arguments{
\item{datalist}{A list used to calculate the (S)IVE-DR estimator including treatment named \code{a}, observed time named \code{obs.t}, censoring indicator (0, censored) named \code{delta}, and baseline covariates used to assign treatment named \code{l}.}
}
\value{
A list including the probability of receiving instrument given baseline covariates named \code{fzl}, the probability of receiving treatment given baseline covariates and instrument equaling 0 named \code{fal0}, the probability of receiving treatment given baseline covariates and instrument equaling 1 named \code{fal1}, and the difference between \code{fal1} and \code{fal0} named \code{deltal}.
}
\description{
Logistic regression for observed treatment and instrument used for the (S)IVE-DR estimator.
}
\details{
More details can be found in references, \code{\link[otrKM]{IVEDR}}, and \code{\link[otrKM]{Genetic.optim.IVEDR}}.
}
\references{
{Xia, J., Zhan, Z., Zhang, J. (2022) An anti-confounding method for estimating optimal regime in a survival context using instrumental variable. Under Review.}
}
