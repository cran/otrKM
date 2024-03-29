#' Given a predetermined t0, estimate the optimal treatment regime by maximizing t0-year survival probability based on the (S)IWKME estimator.
#'
#' @title The optimal treatment regime based on the (S)IWKME estimator.
#' @param datalist A list used to calculate the (S)IWKME estimator including treatment named \code{a}, observed time named \code{obs.t}, censoring indicator (0, censored) named \code{delta}, and baseline covariates used to assign treatment named \code{l}. Notice that all the data in the datalist should be ordered by observed time.
#' @param ps A list including the probability of receiving treatment given baseline covariates named \code{fal}. \code{\link[otrKM]{Fps.IWKME}} can produce \code{ps} by positing logistic model.
#' @param t0 A predetermined time.
#' @param smooth A logic variable indicating wether a smoothed version should be used.
#'
#' @return A numeric vector in which the last number is the estimated optimal t0-year survival probability and the others are the estimated parameter of the optimal treatment regime.
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
#' simulation=simulation[order(simulation$Survival),]
#' 
#' # convert the data into a datalist
#' datalist=list(z=simulation$Instrument,a=simulation$Treatment,
#'               obs.t=simulation$Survival,delta=simulation$Status,
#'               l=cbind(simulation$Covariate1,simulation$Covariate2))
#' 
#' # predetermined t0
#' t0=5
#' 
#' # calculate ps
#' ps=Fps.IWKME(datalist)
#' 
#' Genetic.optim.IWKME(datalist, ps, t0, smooth=TRUE)
Genetic.optim.IWKME <- function(datalist, ps, t0, smooth=TRUE) {
  Neta=ncol(datalist$l)+1
  fn <- IWKME
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