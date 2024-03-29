#' Given a predetermined t0 and eta, calculate t0-year potential survival probability based on the (S)IWKME estimator.
#'
#' @title The (S)IWKME estimator.
#' @param eta The parameters of the regime.
#' @param datalist A list used to calculate the (S)IWKME estimator including treatment named \code{a}, observed time named \code{obs.t}, censoring indicator (0, censored) named \code{delta}, and baseline covariates used to assign treatment named \code{l}. \code{\link[otrKM]{Fps.IWKME}} can produce \code{ps} by positing logistic model. Notice that all the data in the datalist should be ordered by observed time.
#' @param ps A list including the probability of receiving treatment given baseline covariates named \code{fal}.
#' @param t0 A predetermined time.
#' @param smooth A logic variable indicating wether a smoothed estimator should be used.
#'
#' @return Estimated potential survival probability given eta and t0.
#' @export
#' 
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
#' simulation=simulation[order(simulation$Survival),]
#' 
#' # convert the data into a datalist
#' datalist=list(z=simulation$Instrument,a=simulation$Treatment,
#'               obs.t=simulation$Survival,delta=simulation$Status,
#'               l=cbind(simulation$Covariate1,simulation$Covariate2))
#' 
#' # calculate ps
#' ps=Fps.IWKME(datalist)
#' 
#' # predetermined t0 and eta
#' t0=5
#' eta=c(1,2,3)
#' 
#' IWKME(eta, datalist, ps, t0, smooth=TRUE)
IWKME <- function(eta, datalist, ps, t0, smooth=TRUE) {
  # extract variable from datalist
  z <- datalist$z
  a <- datalist$a
  obs.t <- datalist$obs.t
  delta <- datalist$delta
  l <- datalist$l
  fal <- ps$fal

  # regime
  if (!smooth) {
    regime <- as.numeric((cbind(1, l) %*% eta) >= 0)
  } else {
    sd.etal <- sd(c(cbind(1, l) %*% eta))
    if (!is.finite(sd.etal)) {
      return(-1000)
    }
    if (sd.etal > 0) eta <- eta / sd.etal else eta <- c(ifelse(eta[1] >= 0, 1, -1), rep(0, ncol(l)))
    regime <- pnorm(c(cbind(1, l) %*% eta) / ((length(a) / 4)^(-1 / 3)))
  }

  # weight
  w <- (a * regime + (1 - a) * (1 - regime)) / sapply(a * fal + (1 - a) * (1 - fal), clipp)
  w <- w[order(obs.t)]

  # numerator and denominator
  delta.rank <- delta[order(obs.t)]
  num <- w * delta.rank
  den <- sapply(sum(w) - cumsum(w) + w, clipp)

  # survival function
  Nevent <- sum(obs.t <= t0)
  return(prod(1 - num[1:Nevent] / den[1:Nevent]))
}