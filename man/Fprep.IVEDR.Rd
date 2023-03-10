% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{Fprep.IVEDR}
\alias{Fprep.IVEDR}
\title{Cox proportional hazards model for eta-free terms in the (S)IVE-DR estimator.}
\usage{
Fprep.IVEDR(datalist, ps, t0)
}
\arguments{
\item{datalist}{A list used to calculate the (S)IVE-DR estimator including treatment named \code{a}, observed time named \code{obs.t}, censoring indicator (0, censored) named \code{delta}, and baseline covariates used to assign treatment named \code{l}.}

\item{ps}{A list including the probability of receiving instrument given baseline covariates named \code{fzl}, the probability of receiving treatment given baseline covariates and instrument equaling 0 named \code{fal0}, the probability of receiving treatment given baseline covariates and instrument equaling 1 named \code{fal1}, and the difference between fal1 and fal0 named \code{deltal}. \code{\link[otrKM]{Fps.IVEDR}} can produce \code{ps} by positing logistic model.}

\item{t0}{A predetermined t.}
}
\value{
A list including estimates \eqn{\hat{\gamma}_1(\boldsymbol{L};s)}{hat.gamma.1(L;s)} with treatment all to 1 named \code{gamma.num.1} and all to 0 named \code{gamma.num.0}, \eqn{\hat{\gamma}_1'(\boldsymbol{L};s)}{hat.gamma.1'(L;s)} with treatment all to 1 named \code{gammaa.num.1} and all to 0 named \code{gammaa.num.0}, \eqn{\hat{\gamma}_2(\boldsymbol{L};s)}{hat.gamma.2(L;s)} with treatment all to 1 named \code{gamma.den.1} and all to 0 named \code{gamma.den.0}, and \eqn{\hat{\gamma}_2'(\boldsymbol{L};s)}{hat.gamma.2'(L;s)} with treatment all to 1 named \code{gammaa.den.1} and all to 0 named \code{gammaa.den.0}; \code{gamma.num.1} and the others are matrix with ordered observed time as rows and patients as columns. More details can be found in references.
}
\description{
Cox proportional hazards model for eta-free terms in the (S)IVE-DR estimator.
}
\details{
More details can be found in references, \code{\link[otrKM]{IVEDR}}, and \code{\link[otrKM]{Genetic.optim.IVEDR}}.
}
\references{
{Xia, J., Zhan, Z., Zhang, J. (2022) An anti-confounding method for estimating optimal regime in a survival context using instrumental variable. Under Review.}
}
