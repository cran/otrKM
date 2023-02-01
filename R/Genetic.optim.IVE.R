#' Given a predetermined t0, estimate the optimal treatment regime by maximizing t0-year survival probability based on the (S)IVE estimator.
#'
#' @title The optimal treatment regime based on the (S)IVE estimator.
#' @param datalist A list used to calculate the (S)IVE estimator including instrument named z, treatment named \code{a}, observed time named \code{obs.t}, censoring indicator (0, censored) named \code{delta}, and baseline covariates used to assign treatment named \code{l}.
#' @param ps A list including the probability of receiving instrument given baseline covariates named \code{fzl} and the difference between fal1 and fal0 named \code{deltal}, where fal0 denotes the probability of receiving treatment given baseline covariates and instrument equaling 0, and fal1 denotes the probability of receiving treatment given baseline covariates and instrument equaling 1. \code{\link[otrKM]{Fps.IVE}} can produce \code{ps} by positing logistic model.
#' @param t0 A predetermined time.
#' @param smooth A logic variable indicating wether a smoothed version should be used.
#'
#' @return A numeric vector in which the last number is the estimated optimal t0-year survival probability and others are the estimated parameter of the optimal treatment regime.
#' @export
#' 
#' @details More details can be found in references.
#' @references 
#' {Xia, J., Zhan, Z., Zhang, J. (2022) An anti-confounding method for estimating optimal regime in a survival context using instrumental variable. Under Review.}
#' 
#' @import survival
#' @import rgenoud
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
#' # predetermined t0
#' t0=5
#' 
#' Genetic.optim.IVE(datalist, ps, t0, smooth=TRUE)
Genetic.optim.IVE <- function(datalist, ps, t0, smooth=TRUE) {
  Neta=ncol(datalist$l)+1
  fn <- IVE
  temp <- genoud(
    fn = fn, datalist = datalist, ps = ps, t0 = t0, smooth = smooth,
    nvars = Neta,
    Domains = cbind(rep(-1, Neta), rep(1, Neta)),
    starting.values = rep(0.001, Neta),
    max = TRUE,
    pop.size = 20,
    wait.generations = 10,
    print.level = 0,
    unif.seed = 1107,
    int.seed = 0130
  )

  eta.est <- temp$par
  etahat <- eta.est / sqrt(sum(eta.est^2))
  Shat.etahat <- temp$value

  return(c(etahat, Shat.etahat))
}
