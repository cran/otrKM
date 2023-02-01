#' Given a predetermined t0, estimate the optimal treatment regime by maximizing t0-year survival probability based on the (S)AIPS estimator.
#'
#' @title The optimal treatment regime based on the (S)AIPS estimator.
#' @param datalist A list used to calculate the (S)AIPS estimator including treatment named \code{a}, observed time named \code{obs.t}, censoring indicator (0, censored) named \code{delta}, and baseline covariates used to assign treatment named \code{l}.
#' @param ps A list including the probability of receiving treatment given baseline covariates named \code{fal}. \code{\link[otrKM]{Fps.AIPS}} can produce \code{ps} by positing logistic model.
#' @param prep A list including the augmented terms in the numerator with treatment all to 1 named \code{gamma.num.1} and all to 0 named \code{gamma.num.0} and in the denominator with treatment all to 1 named \code{gamma.den.1} and all to 0 named \code{gamma.den.0}; \code{gamma.num.1} and the others are matrix with ordered observed time as rows and patients as columns. \code{\link[otrKM]{Fprep.AIPS}} can produce \code{prep} by positing Cox proportional hazards model.
#' @param t0 A predetermined time.
#' @param smooth A logic variable indicating wether a smoothed version should be used.
#'
#' @return A numeric vector in which the last number is the estimated optimal t0-year survival probability and others are the estimated parameters of the optimal treatment regime.
#' @export
#' 
#' @details More details can be found in references.
#' @references 
#' {Jiang, R., Lu, W., Song, R., and Davidian, M. (2017) On estimation of optimal treatment regimes for maximizing t‐year survival probability. \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, \bold{79:} 1165-1185. DOI:10.1111/rssb.12201}
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
#' # predetermined t0
#' t0=5
#' 
#' # calculate ps and prep
#' ps=Fps.AIPS(datalist)
#' prep=Fprep.AIPS(datalist, t0)
#' 
#' Genetic.optim.AIPS(datalist, ps, prep, t0, smooth=TRUE)
Genetic.optim.AIPS <- function(datalist, ps, prep, t0, smooth=TRUE) {
  Neta=ncol(datalist$l)+1
  fn <- AIPS
  temp <- genoud(
    fn = fn, datalist = datalist, ps = ps, prep = prep, t0 = t0, smooth = smooth,
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