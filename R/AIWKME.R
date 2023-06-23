#' Given a predetermined t0 and eta, calculate t0-year potential survival probability based on the (S)AIWKME estimator.
#' 
#' @title The (S)AIWKME estimator.
#' @param eta The parameters of the regime.
#' @param datalist A list used to calculate the (S)AIWKME estimator including treatment named \code{a}, observed time named \code{obs.t}, censoring indicator (0, censored) named \code{delta}, and baseline covariates used to assign treatment named \code{l}.
#' @param ps A list including the probability of receiving treatment given baseline covariates named \code{fal}. \code{\link[otrKM]{Fps.AIWKME}} can produce \code{ps} by positing logistic model.
#' @param prep A list including the augmented terms in the numerator with treatment all to 1 named \code{gamma.num.1} and all to 0 named \code{gamma.num.0} and in the denominator with treatment all to 1 named \code{gamma.den.1} and all to 0 \code{named gamma.den.0}; \code{gamma.num.1} and the others are matrix with ordered observed time as rows and patients as columns. \code{\link[otrKM]{Fprep.AIWKME}} can produce \code{prep} by positing Cox proportional hazards model.
#' @param t0 A predetermined time. 
#' @param smooth A logic variable indicating wether a smoothed estimator should be used.
#'
#' @return Estimated potential survival probability given eta and t0.
#' @export
#' @details More details can be found in references.
#' @references 
#' {Jiang, R., Lu, W., Song, R., and Davidian, M. (2017) On estimation of optimal treatment regimes for maximizing t‐year survival probability. \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, \bold{79:} 1165-1185. DOI:10.1111/rssb.12201}
#' 
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
#' # predetermined t0 and eta
#' t0=5
#' eta=c(1,2,3)
#' 
#' # calculate ps and prep
#' ps=Fps.AIWKME(datalist)
#' prep=Fprep.AIWKME(datalist, t0)
#' 
#' AIWKME(eta, datalist, ps, prep, t0, smooth=TRUE)
AIWKME <- function(eta, datalist, ps, prep, t0, smooth=TRUE) {
  # extract variable from datalist
  a <- datalist$a
  obs.t <- datalist$obs.t
  delta <- datalist$delta
  l <- datalist$l
  fal <- ps$fal
  gamma.num.1 <- prep$gamma.num.1
  gamma.num.0 <- prep$gamma.num.0
  gamma.den.1 <- prep$gamma.den.1
  gamma.den.0 <- prep$gamma.den.0

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
  w0 <- (1 - a) / sapply(fal*a+(1 - fal)*(1-a), clipp)
  w1 <- a / sapply(fal*a+(1-fal)*(1-a), clipp)

  delta.rank <- delta[order(obs.t)]
  w <- (w1 * regime + w0 * (1 - regime))[order(obs.t)]
  Nevent <- sum(obs.t <= t0)

  # numerator and denominator
  num <- (w * delta.rank)[1:Nevent] + gamma.num.1 %*% (regime * (1 - w1)) + gamma.num.0 %*% ((1 - regime) * (1 - w0))

  den <- (sum(w) - cumsum(w) + w)[1:Nevent] + gamma.den.1 %*% (regime * (1 - w1)) + gamma.den.0 %*% ((1 - regime) * (1 - w0))

  # survival probability
  return(prod(1 - num / den))
}