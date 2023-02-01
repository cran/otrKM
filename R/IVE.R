#' Given a predetermined t0 and eta, calculate t0-year potential survival probability based on the (S)IVE estimator.
#'
#' @title The (S)IVE estimator.
#' @param eta The parameters of the regime.
#' @param datalist A list used to calculate the (S)IVE estimator including treatment named \code{a}, observed time named \code{obs.t}, censoring indicator (0, censored) named \code{delta}, and baseline covariates used to assign treatment named \code{l}.
#' @param ps A list including the probability of receiving instrument given baseline covariates named \code{fzl} and the difference between fal1 and fal0 named \code{deltal}, where fal0 denotes the probability of receiving treatment given baseline covariates and instrument equaling 0, and fal1 denotes the probability of receiving treatment given baseline covariates and instrument equaling 1. \code{\link[otrKM]{Fps.IVE}} can produce \code{ps} by positing logistic model.
#' @param t0 A predetermined time.
#' @param smooth A logic variable indicating wether a smoothed version should be used.
#'
#' @return Estimated potential survival probability given eta and t0.
#' @export
#' @details More details can be found in references.
#' @references 
#' {Xia, J., Zhan, Z., Zhang, J. (2022) An anti-confounding method for estimating optimal regime in a survival context using instrumental variable. Under Review.}
#' 
#' @import survival
#' @import rgenoud
#' @import stats
#'
#' @examples
#' # load data
#' data(simulation)
#' 
#' # convert the data into a datalist
#' datalist=list(z=simulation$Instrument,a=simulation$Treatment,
#'               obs.t=simulation$Survival,delta=simulation$Status,
#'               l=cbind(simulation$Covariate1,simulation$Covariate2))
#' 
#' # calculate ps and prep
#' ps=Fps.IVE(datalist)
#' 
#' # predetermined eta
#' t0=5
#' eta=c(1,2,3)
#' 
#' IVE(eta, datalist, ps, t0, smooth=TRUE)
IVE <- function(eta, datalist, ps, t0, smooth=TRUE) {
  # extract variable from datalist
  z <- datalist$z
  a <- datalist$a
  obs.t <- datalist$obs.t
  delta <- datalist$delta
  l <- datalist$l
  fzl <- ps$fzl
  deltal <- ps$deltal

  # regime
  if (!smooth) {
    regime <- as.numeric((cbind(1, l) %*% eta) >= 0)
  } else {
    sd.etal <- sd(c(cbind(1, l) %*% eta))
    if (!is.finite(sd.etal)) return(-1000)
    if (sd.etal > 0) eta <- eta / sd.etal else eta <- c(ifelse(eta[1] >= 0, 1, -1), rep(0, ncol(l)))
    regime <- pnorm(c(cbind(1, l) %*% eta) / ((length(a) / 4)^(-1 / 3)))
  }

  # weight
  w <- (2 * a - 1) * (2 * z - 1) * (a * regime + (1 - a) * (1 - regime)) / sapply((z * fzl + (1 - z) * (1 - fzl)) * deltal, clipp)
  w <- w[order(obs.t)]

  # numerator and denominator
  delta.rank <- delta[order(obs.t)]
  num <- w * delta.rank
  den <- sum(w) - cumsum(w) + w

  # survival function
  Nevent <- sum(obs.t <= t0)
  surv <- prod(1 - num[1:Nevent] / den[1:Nevent])
  return(surv)
}
