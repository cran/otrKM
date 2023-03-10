% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{Fps.AIPS}
\alias{Fps.AIPS}
\title{Logistic regression for observed treatment used for the (S)AIPS estimator.}
\usage{
Fps.AIPS(datalist)
}
\arguments{
\item{datalist}{A list used to calculate the (S)AIPS estimator including treatment named \code{a}, observed time named \code{obs.t}, censoring indicator (0, censored) named \code{delta}, and baseline covariates used to assign treatment named \code{l}.}
}
\value{
A list including the probability of receiving treatment given baseline covariates named \code{fal}.
}
\description{
Logistic regression for observed treatment used for the (S)AIPS estimator.
}
\details{
More details can be found in references, \code{\link[otrKM]{AIPS}}, and \code{\link[otrKM]{Genetic.optim.AIPS}}.
}
\references{
{Jiang, R., Lu, W., Song, R., and Davidian, M. (2017) On estimation of optimal treatment regimes for maximizing t‐year survival probability. \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, \bold{79:} 1165-1185. DOI:10.1111/rssb.12201}
}
