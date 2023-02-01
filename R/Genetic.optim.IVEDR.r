#' Given a predetermined t0, estimate the optimal treatment regime by maximizing t0-year survival probability based on the (S)IVE-DR estimator.
#' 
#' @title The optimal treatment regime based on the (S)IVE-DR estimator.
#' @param datalist A list used to calculate the (S)IVE-DR estimator including instrument named \code{z}, treatment named \code{a}, observed time named \code{obs.t}, censoring indicator (0, censored) named \code{delta}, and baseline covariates used to assign treatment named \code{l}.
#' @param ps A list including the probability of receiving instrument given baseline covariates named \code{fzl}, the probability of receiving treatment given baseline covariates and instrument equaling 0 named \code{fal0}, the probability of receiving treatment given baseline covariates and instrument equaling 1 named \code{fal1}, and the difference between fal1 and fal0 named \code{deltal}. \code{\link[otrKM]{Fps.IVEDR}} can produce \code{ps} by positing logistic model.
#' @param prep A list including estimates \eqn{\hat{\gamma}_1(\boldsymbol{L};s)}{hat.gamma.1(L;s)} with treatment all to 1 named \code{gamma.num.1} and all to 0 named \code{gamma.num.0}, \eqn{\hat{\gamma}_1'(\boldsymbol{L};s)}{hat.gamma.1'(L;s)} with treatment all to 1 named \code{gammaa.num.1} and all to 0 named \code{gammaa.num.0}, \eqn{\hat{\gamma}_2(\boldsymbol{L};s)}{hat.gamma.2(L;s)} with treatment all to 1 named \code{gamma.den.1} and all to 0 named \code{gamma.den.0}, and \eqn{\hat{\gamma}_2'(\boldsymbol{L};s)}{hat.gamma.2'(L;s)} with treatment all to 1 named \code{gammaa.den.1} and all to 0 named \code{gammaa.den.0}; \code{gamma.num.1} and the others are matrix with ordered observed time as rows and patients as columns. More details can be found in references. \code{\link[otrKM]{Fprep.IVEDR}} can produce \code{prep} by positing Cox proportional hazards model.
#' @param t0 A predetermined time to point out that t0-year survival probability is our estimate
#' @param smooth A logic variable indicating wether a smoothed version should be used.
#'
#' @return A numeric vector in which the last number is the estimated optimal t0-year survival probability and the others are the estimated parameter of the optimal treatment regime.
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
#' #' # predetermined t0
#' t0=1
#' 
#' # calculate ps and prep
#' ps=Fps.IVEDR(datalist)
#' prep=Fprep.IVEDR(datalist, ps, t0)
#' 
#' 
#' Genetic.optim.IVEDR(datalist, ps, prep, t0, smooth=TRUE)
Genetic.optim.IVEDR <- function(datalist, ps, prep, t0, smooth=TRUE){
  Neta=ncol(datalist$l)+1
  fn <- IVEDR
  temp <- genoud(fn=fn, datalist=datalist, ps=ps, t0=t0, smooth=smooth, prep=prep,
                 nvars=Neta,
                 Domains=cbind(rep(-1,Neta),rep(1,Neta)),
                 starting.values=rep(0.001,Neta),

                 max=TRUE,
                 pop.size=20,
                 wait.generations=10,
                 print.level=0,

                 unif.seed=1107,
                 int.seed=0130
  )

  eta.est <- temp$par
  etahat <- eta.est/sqrt(sum(eta.est^2))
  Shat.etahat <- temp$value

  return(c(etahat, Shat.etahat))
}