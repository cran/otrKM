#' Limit the number not to be too large or too small.
#' @title clip function.
#' @param x A vector or matrix.
#' @export
#' @return A vector or matrix same as the input.
clipp <- function(x) {
  return(sign(x) * min(max(abs(x), 0.00001), 100000))
}

#' Logistic regression for observed treatment used for the (S)IPS estimator.
#' @param datalist A list used to calculate the (S)IPS estimator including treatment named \code{a}, observed time named \code{obs.t}, censoring indicator (0, censored) named \code{delta}, and baseline covariates used to assign treatment named \code{l}.
#' @import stats
#' @return A list including the probability of receiving treatment given baseline covariates named \code{fal}.
#' @export
#' @details More details can be found in references, \code{\link[otrKM]{IPS}}, and \code{\link[otrKM]{Genetic.optim.IPS}}.
#' @references 
#' {Jiang, R., Lu, W., Song, R., and Davidian, M. (2017) On estimation of optimal treatment regimes for maximizing t‐year survival probability. \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, \bold{79:} 1165-1185. DOI:10.1111/rssb.12201}
Fps.IPS <- function(datalist) {
  l <- datalist$l
  a <- datalist$a

  f.a.l <- glm(a ~ l, family = binomial(link = "logit"))
  fal <- predict.glm(f.a.l, type = "response")

  return(list(fal = fal))
}

#' Logistic regression for observed treatment used for the (S)AIPS estimator.
#' @param datalist A list used to calculate the (S)AIPS estimator including treatment named \code{a}, observed time named \code{obs.t}, censoring indicator (0, censored) named \code{delta}, and baseline covariates used to assign treatment named \code{l}.
#' @import stats
#' @return A list including the probability of receiving treatment given baseline covariates named \code{fal}.
#' @export
#' @details More details can be found in references, \code{\link[otrKM]{AIPS}}, and \code{\link[otrKM]{Genetic.optim.AIPS}}.
#' @references 
#' {Jiang, R., Lu, W., Song, R., and Davidian, M. (2017) On estimation of optimal treatment regimes for maximizing t‐year survival probability. \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, \bold{79:} 1165-1185. DOI:10.1111/rssb.12201}
Fps.AIPS <- function(datalist) {
  l <- datalist$l
  a <- datalist$a

  f.a.l <- glm(a ~ l, family = binomial(link = "logit"))
  fal <- predict.glm(f.a.l, type = "response")

  return(list(fal = fal))
}

#' Cox proportional hazards model for eta-free terms in the (S)AIPS estimator.
#' @param datalist A list used to calculate the (S)AIPS estimator including treatment named \code{a}, observed time named \code{obs.t}, censoring indicator (0, censored) named \code{delta}, and baseline covariates used to assign treatment named \code{l}.
#' @param t0 A predetermined t.
#' @import stats
#' @import survival
#' @return A list including the augmented terms in the numerator with treatment all to 1 named \code{gamma.num.1} and all to 0 named \code{gamma.num.0} and in the denominator with treatment all to 1 named \code{gamma.den.1} and all to 0 named \code{gamma.den.0}; \code{gamma.num.1} and the others are matrix with ordered observed time as rows and patients as columns. More details can be found in references.
#' @export
#' @details More details can be found in references, \code{\link[otrKM]{AIPS}}, and \code{\link[otrKM]{Genetic.optim.AIPS}}.
#' @references 
#' {Jiang, R., Lu, W., Song, R., and Davidian, M. (2017) On estimation of optimal treatment regimes for maximizing t‐year survival probability. \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, \bold{79:} 1165-1185. DOI:10.1111/rssb.12201}
Fprep.AIPS <- function(datalist, t0) {
  # extract variable from datalist
  z <- datalist$z
  a <- datalist$a
  obs.t <- datalist$obs.t
  delta <- datalist$delta
  l <- datalist$l
  Nevent <- sum(obs.t <= t0)
  rank <- order(obs.t)

  # estimate the survival function of the censoring time
  fit.C <- survfit(Surv(obs.t, 1 - delta) ~ 1)
  sur.C <- rep(fit.C$surv, times = (fit.C$n.event + fit.C$n.censor))[1:Nevent]

  # estimate the hazards function
  cox <- coxph(Surv(obs.t, delta) ~ l * a, ties = "breslow")

  # estimate baseline hazard function
  cox.a <- c(exp(cbind(l, a, l * a) %*% cox$coefficients))
  cox.a <- cox.a[rank]
  delta.rank <- delta[rank]
  den.cox <- sum(cox.a) - cumsum(cox.a) + cox.a
  Lambda0 <- cumsum(delta.rank / den.cox)[1:Nevent] # baseline cumulative hazard function
  bhest <- c(Lambda0[1], diff(Lambda0))[1:Nevent] # baseline hazard function

  cox.predict.0 <- c(exp(cbind(l, 0, l * 0) %*% cox$coefficients))
  cox.predict.1 <- c(exp(cbind(l, 1, l) %*% cox$coefficients))

  surv.T.0 <- exp(-Lambda0 %o% cox.predict.0)
  surv.T.1 <- exp(-Lambda0 %o% cox.predict.1)


  gamma.num.1 <- sur.C * surv.T.1 * (bhest %o% cox.predict.1)
  gamma.num.0 <- sur.C * surv.T.0 * (bhest %o% cox.predict.0)

  gamma.den.1 <- sur.C * surv.T.1
  gamma.den.0 <- sur.C * surv.T.0

  return(list(
    gamma.num.1 = gamma.num.1,
    gamma.num.0 = gamma.num.0,
    gamma.den.1 = gamma.den.1,
    gamma.den.0 = gamma.den.0
  ))
}

#' Logistic regression for observed treatment and instrument used for the (S)IVE estimator. 
#' @param datalist A list used to calculate the (S)IVE estimator including treatment named \code{a}, observed time named \code{obs.t}, censoring indicator (0, censored) named \code{delta}, and baseline covariates used to assign treatment named \code{l}.
#' @import stats
#' @return A list including the probability of receiving instrument given baseline covariates named \code{fzl} and the difference between fal1 and fal0 named \code{deltal}, where fal0 denotes the probability of receiving treatment given baseline covariates and instrument equaling 0, and fal1 denotes the probability of receiving treatment given baseline covariates and instrument equaling 1.
#' @details More details can be found in references, \code{\link[otrKM]{IVE}}, and \code{\link[otrKM]{Genetic.optim.IVE}}.
#' @references 
#' {Xia, J., Zhan, Z., Zhang, J. (2022) An anti-confounding method for estimating optimal regime in a survival context using instrumental variable. Under Review.}
#' @export
Fps.IVE <- function(datalist) {
  l <- datalist$l
  z <- datalist$z
  a <- datalist$a

  f.z.l <- glm(z ~ l, family = binomial(link = "logit"))
  fzl <- predict.glm(f.z.l, type = "response")

  f.a.lz <- glm(a ~ l + z, family = binomial(link = "logit"))
  fal1 <- predict.glm(f.a.lz, newdata = data.frame(l, z = 1), type = "response")
  fal0 <- predict.glm(f.a.lz, newdata = data.frame(l, z = 0), type = "response")
  deltal <- fal1 - fal0

  return(list(fzl = fzl, deltal = deltal, fal0 = fal0, fal1=fal1))
}

#' Logistic regression for observed treatment and instrument used for the (S)IVE-DR estimator.
#' @param datalist A list used to calculate the (S)IVE-DR estimator including treatment named \code{a}, observed time named \code{obs.t}, censoring indicator (0, censored) named \code{delta}, and baseline covariates used to assign treatment named \code{l}.
#' @import stats
#' @return A list including the probability of receiving instrument given baseline covariates named \code{fzl}, the probability of receiving treatment given baseline covariates and instrument equaling 0 named \code{fal0}, the probability of receiving treatment given baseline covariates and instrument equaling 1 named \code{fal1}, and the difference between \code{fal1} and \code{fal0} named \code{deltal}.
#' @export
#' @details More details can be found in references, \code{\link[otrKM]{IVEDR}}, and \code{\link[otrKM]{Genetic.optim.IVEDR}}.
#' @references 
#' {Xia, J., Zhan, Z., Zhang, J. (2022) An anti-confounding method for estimating optimal regime in a survival context using instrumental variable. Under Review.}
Fps.IVEDR <- function(datalist){
  l=datalist$l
  z=datalist$z
  a=datalist$a

  f.z.l=glm(z~l,family = binomial(link='logit'))
  fzl=predict.glm(f.z.l,type='response')

  f.a.lz=glm(a~l+z,family = binomial(link='logit'))
  fal1=predict.glm(f.a.lz,newdata=data.frame(l,z=1),type='response')
  fal0=predict.glm(f.a.lz,newdata=data.frame(l,z=0),type='response')
  deltal=fal1-fal0

  return(list(fzl=fzl, deltal=deltal,fal0=fal0, fal1=fal1))
}

#' Cox proportional hazards model for eta-free terms in the (S)IVE-DR estimator.
#' @param datalist A list used to calculate the (S)IVE-DR estimator including treatment named \code{a}, observed time named \code{obs.t}, censoring indicator (0, censored) named \code{delta}, and baseline covariates used to assign treatment named \code{l}.
#' @param ps A list including the probability of receiving instrument given baseline covariates named \code{fzl}, the probability of receiving treatment given baseline covariates and instrument equaling 0 named \code{fal0}, the probability of receiving treatment given baseline covariates and instrument equaling 1 named \code{fal1}, and the difference between fal1 and fal0 named \code{deltal}. \code{\link[otrKM]{Fps.IVEDR}} can produce \code{ps} by positing logistic model.
#' @param t0 A predetermined t.
#' 
#' @import stats
#' @import survival
#' @return A list including estimates \eqn{\hat{\gamma}_1(\boldsymbol{L};s)}{hat.gamma.1(L;s)} with treatment all to 1 named \code{gamma.num.1} and all to 0 named \code{gamma.num.0}, \eqn{\hat{\gamma}_1'(\boldsymbol{L};s)}{hat.gamma.1'(L;s)} with treatment all to 1 named \code{gammaa.num.1} and all to 0 named \code{gammaa.num.0}, \eqn{\hat{\gamma}_2(\boldsymbol{L};s)}{hat.gamma.2(L;s)} with treatment all to 1 named \code{gamma.den.1} and all to 0 named \code{gamma.den.0}, and \eqn{\hat{\gamma}_2'(\boldsymbol{L};s)}{hat.gamma.2'(L;s)} with treatment all to 1 named \code{gammaa.den.1} and all to 0 named \code{gammaa.den.0}; \code{gamma.num.1} and the others are matrix with ordered observed time as rows and patients as columns. More details can be found in references.
#' 
#' 
#' @export
#' @details More details can be found in references, \code{\link[otrKM]{IVEDR}}, and \code{\link[otrKM]{Genetic.optim.IVEDR}}.
#' @references 
#' {Xia, J., Zhan, Z., Zhang, J. (2022) An anti-confounding method for estimating optimal regime in a survival context using instrumental variable. Under Review.}
Fprep.IVEDR <- function(datalist, ps, t0){
  # extract variable from list
  z=datalist$z
  a=datalist$a
  obs.t=datalist$obs.t
  delta=datalist$delta
  l=datalist$l
  fal0=ps$fal0
  fal1=ps$fal1
  deltal=ps$deltal

  Nevent=sum(obs.t<=t0)
  rank=order(obs.t)

  # estimate the survival function of the censoring time
  fit.C <- survfit(Surv(obs.t,1-delta)~1)
  sur.C=rep(fit.C$surv, times = (fit.C$n.event + fit.C$n.censor))[1:Nevent]

  # estimate the hazards function
  cox <- coxph(Surv(obs.t, delta)~l*a+z,ties='breslow')

  # estimate baseline hazard function
  cox.a=c(exp(cbind(l,a,z,l*a)%*%cox$coefficients))
  cox.a=cox.a[rank]
  delta.rank=delta[rank]
  den.cox <- sum(cox.a) - cumsum(cox.a) + cox.a
  Lambda0 <- cumsum(delta.rank / den.cox)[1:Nevent] # baseline cumulative hazard function
  bhest <- c(Lambda0[1], diff(Lambda0))[1:Nevent] # baseline hazard function

  # In the following, the matrix will be ordered observed t times patients
  cox.predict.11<-c(exp(cbind(l,1,1,l)%*%cox$coefficients))
  cox.predict.10<-c(exp(cbind(l,1,0,l)%*%cox$coefficients))

  # the first 1 represent a=1, while the second 1 represent z=1
  surv.T.11=exp(-Lambda0%o%cox.predict.11)
  surv.T.10=exp(-Lambda0%o%cox.predict.10)

  gammaa.num.1=t(t(sur.C*surv.T.10*(bhest%o%cox.predict.10))*fal0)

  gamma.num.1=t(t(sur.C*surv.T.11*(bhest%o%cox.predict.11))*fal1/sapply(deltal,clipp))-t(t(sur.C*surv.T.10*(bhest%o%cox.predict.10))*fal0/sapply(deltal,clipp))

  gammaa.den.1=t(t(sur.C*surv.T.10)*fal0)

  gamma.den.1=t(t(sur.C*surv.T.11)*fal1/sapply(deltal,clipp))-t(t(sur.C*surv.T.10)*fal0/sapply(deltal,clipp))

  cox.predict.01<-c(exp(cbind(l,0,1,l*0)%*%cox$coefficients))
  cox.predict.00<-c(exp(cbind(l,0,0,l*0)%*%cox$coefficients))

  surv.T.01=exp(-Lambda0%o%cox.predict.01)
  surv.T.00=exp(-Lambda0%o%cox.predict.00)

  gammaa.num.0=-t(t(sur.C*surv.T.00*(bhest%o%cox.predict.00))*(1-fal0))

  gamma.num.0=-t(t(sur.C*surv.T.01*(bhest%o%cox.predict.01))*(1-fal1)/sapply(deltal,clipp))+t(t(sur.C*surv.T.00*(bhest%o%cox.predict.00))*(1-fal0)/sapply(deltal,clipp))

  gammaa.den.0=-t(t(sur.C*surv.T.00)*(1-fal0))

  gamma.den.0=-t(t(sur.C*surv.T.01)*(1-fal1)/sapply(deltal,clipp))+t(t(sur.C*surv.T.00)*(1-fal0)/sapply(deltal,clipp))

  return(list(  gammaa.num.1=gammaa.num.1,
                gammaa.num.0=gammaa.num.0,
                gamma.num.1=gamma.num.1,
                gamma.num.0=gamma.num.0,
                gammaa.den.1=gammaa.den.1,
                gammaa.den.0=gammaa.den.0,
                gamma.den.1=gamma.den.1,
                gamma.den.0=gamma.den.0))
}