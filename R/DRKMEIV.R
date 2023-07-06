#' Given a predetermined t0 and eta, calculate t0-year potential survival probability based on the (S)DRKMEIV estimator.
#' 
#' @title The (S)DRKMEIV estimator.
#' @param eta The parameters of the regime.
#' @param datalist A list used to calculate the (S)DRKMEIV estimator including treatment named \code{a}, observed time named \code{obs.t}, censoring indicator (0, censored) named \code{delta}, and baseline covariates used to assign treatment named \code{l}. Notice that all the data in the datalist should be ordered by observed time.
#' @param ps A list including the probability of receiving instrument given baseline covariates named \code{fzl}, the probability of receiving treatment given baseline covariates and instrument equaling 0 named \code{fal0}, the probability of receiving treatment given baseline covariates and instrument equaling 1 named \code{fal1}, and the difference between fal1 and fal0 named \code{deltal}. \code{\link[otrKM]{Fps.DRKMEIV}} can produce \code{ps} by positing logistic model.
#' @param prep A list including estimates \eqn{\hat{\gamma}_1(\boldsymbol{L};s)}{hat.gamma.1(L;s)} with treatment all to 1 named \code{gamma.num.1} and all to 0 named \code{gamma.num.0}, \eqn{\hat{\gamma}_1'(\boldsymbol{L};s)}{hat.gamma.1'(L;s)} with treatment all to 1 named \code{gammaa.num.1} and all to 0 named \code{gammaa.num.0}, \eqn{\hat{\gamma}_2(\boldsymbol{L};s)}{hat.gamma.2(L;s)} with treatment all to 1 named \code{gamma.den.1} and all to 0 named \code{gamma.den.0}, and \eqn{\hat{\gamma}_2'(\boldsymbol{L};s)}{hat.gamma.2'(L;s)} with treatment all to 1 named \code{gammaa.den.1} and all to 0 named \code{gammaa.den.0}; \code{gamma.num.1} and the others are matrix with ordered observed time as rows and patients as columns. There are also estimates for the last term of the (S)DRKMEIV estimator. More details can be found in references. \code{\link[otrKM]{Fprep.DRKMEIV}} can produce \code{prep} by positing Cox proportional hazards model.
#' @param t0 A predetermined time.
#' @param smooth A logic variable indicating wether a smoothed version should be used.
#'
#' @return Estimated potential survival probability given eta and t0.
#' @export
#' 
#' @details More details can be found in references.
#' @references 
#' {Xia, J., Zhan, Z., Zhang, J. (2022) Estimating optimal treatment regime in survival contexts using an instrumental variable. Under Review.}
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
#' # predetermined t0 and eta
#' t0=5
#' eta=c(1,2,3)
#' 
#' # calculate ps and prep
#' ps=Fps.DRKMEIV(datalist,t0)
#' prep=Fprep.DRKMEIV(datalist, ps, t0)
#' 
#' DRKMEIV(eta, datalist, ps, prep, t0, smooth=TRUE)
DRKMEIV<-function(eta, datalist, ps, prep, t0, smooth=TRUE){
  # extract variable from datalist
  z=datalist$z
  a=datalist$a
  obs.t=datalist$obs.t
  delta=datalist$delta
  l=datalist$l

  fzl=ps$fzl
  deltal=ps$deltal
  fal0=ps$fal0
  surv.C.1=ps$surv.C.1
  surv.C.0=ps$surv.C.0
  
  gammaa.num.1=prep$gammaa.num.1
  gammaa.num.0=prep$gammaa.num.0
  gamma.num.1=prep$gamma.num.1
  gamma.num.0=prep$gamma.num.0
  gammaa.den.1=prep$gammaa.den.1
  gammaa.den.0=prep$gammaa.den.0
  gamma.den.1=prep$gamma.den.1
  gamma.den.0=prep$gamma.den.0
  mat.num.1=prep$mat.num.1
  mat.num.0=prep$mat.num.0
  mat.den.1=prep$mat.den.1
  mat.den.0=prep$mat.den.0

  # regime
  if (!smooth){
    regime <- as.numeric((cbind(1,l) %*% eta) >= 0)
  }
  else{
    sd.etal <- sd(c(cbind(1,l) %*% eta))
    if (!is.finite(sd.etal)) return(-1000)
    if (sd.etal > 0) eta <- eta/sd.etal else eta <- c(ifelse(eta[1] >= 0, 1, -1), rep(0, ncol(l)))
    regime <- pnorm(c(cbind(1,l) %*% eta)/((length(a)/4)^(-1/3)))
  }

  # weight
  Nevent <- sum(obs.t <= t0)
  delta.rank <- delta[order(obs.t)][1:Nevent]
  dN <- cbind(diag(delta.rank), matrix(0, nrow = Nevent, ncol = length(a)-Nevent))
  Y<-dN
  Y[upper.tri(Y,diag=TRUE)] <- 1

  w1 <- (2 * a - 1) * (2 * z - 1) * a*regime / sapply((z * fzl + (1 - z) * (1 - fzl)) * deltal, clipp)
  w0<- (2 * a - 1) * (2 * z - 1) * (1 - a)*(1-regime) / sapply((z * fzl + (1 - z) * (1 - fzl)) * deltal, clipp)


  # numerator and denominator
  num<-c((dN/surv.C.1)%*%w1+(dN/surv.C.0)%*%w0)-gammaa.num.1%*%(regime*(2*z-1)/sapply((z*fzl+(1-z)*(1-fzl))*deltal,clipp))-gamma.num.1%*%(regime*(a-fal0)*(2*z-1)/sapply((z*fzl+(1-z)*(1-fzl))*deltal,clipp))+rowSums(gamma.num.1%*%regime)-gammaa.num.0%*%((1-regime)*(2*z-1)/sapply((z*fzl+(1-z)*(1-fzl))*deltal,clipp))-gamma.num.0%*%((1-regime)*(a-fal0)*(2*z-1)/sapply((z*fzl+(1-z)*(1-fzl))*deltal,clipp))+rowSums(gamma.num.0%*%(1-regime)) +c(mat.num.1%*%w1+mat.num.0%*%w0)

  den<-c((Y/surv.C.1)%*%w1+(Y/surv.C.0)%*%w0)-gammaa.den.1%*%(regime*(2*z-1)/sapply((z*fzl+(1-z)*(1-fzl))*deltal,clipp))-gamma.den.1%*%(regime*(a-fal0)*(2*z-1)/sapply((z*fzl+(1-z)*(1-fzl))*deltal,clipp))+rowSums(gamma.den.1%*%regime)-gammaa.den.0%*%((1-regime)*(2*z-1)/sapply((z*fzl+(1-z)*(1-fzl))*deltal,clipp))-gamma.den.0%*%((1-regime)*(a-fal0)*(2*z-1)/sapply((z*fzl+(1-z)*(1-fzl))*deltal,clipp))+rowSums(gamma.den.0%*%(1-regime))+c(mat.den.1%*%w1+mat.den.0%*%w0)

  # survival function
  return(prod(1-num/den))
}