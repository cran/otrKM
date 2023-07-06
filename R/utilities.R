#' Limit the number not to be too large or too small.
#' @title clip function.
#' @param x A vector or matrix.
#' @export
#' @return A vector or matrix same as the input.
clipp <- function(x) {
  return(sign(x) * min(max(abs(x), 0.00001), 100000))
}

#' Logistic regression for observed treatment used for the (S)IWKME estimator.
#' @param datalist A list used to calculate the (S)IWKME estimator including treatment named \code{a}, observed time named \code{obs.t}, censoring indicator (0, censored) named \code{delta}, and baseline covariates used to assign treatment named \code{l}. Notice that all the data in the datalist should be ordered by observed time.
#' @import stats
#' @return A list including the probability of receiving treatment given baseline covariates named \code{fal}.
#' @export
#' @details More details can be found in references, \code{\link[otrKM]{IWKME}}, and \code{\link[otrKM]{Genetic.optim.IWKME}}.
#' @references 
#' {Jiang, R., Lu, W., Song, R., and Davidian, M. (2017) On estimation of optimal treatment regimes for maximizing t‐year survival probability. \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, \bold{79:} 1165-1185. DOI:10.1111/rssb.12201}
Fps.IWKME <- function(datalist) {
  l <- datalist$l
  a <- datalist$a

  f.a.l <- glm(a ~ l, family = binomial(link = "logit"))
  fal <- predict.glm(f.a.l, type = "response")

  return(list(fal = fal))
}

#' Logistic regression for observed treatment used for the (S)AIWKME estimator.
#' @param datalist A list used to calculate the (S)AIWKME estimator including treatment named \code{a}, observed time named \code{obs.t}, censoring indicator (0, censored) named \code{delta}, and baseline covariates used to assign treatment named \code{l}. Notice that all the data in the datalist should be ordered by observed time.
#' @import stats
#' @return A list including the probability of receiving treatment given baseline covariates named \code{fal}.
#' @export
#' @details More details can be found in references, \code{\link[otrKM]{AIWKME}}, and \code{\link[otrKM]{Genetic.optim.AIWKME}}.
#' @references 
#' {Jiang, R., Lu, W., Song, R., and Davidian, M. (2017) On estimation of optimal treatment regimes for maximizing t‐year survival probability. \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, \bold{79:} 1165-1185. DOI:10.1111/rssb.12201}
Fps.AIWKME <- function(datalist) {
  l <- datalist$l
  a <- datalist$a

  f.a.l <- glm(a ~ l, family = binomial(link = "logit"))
  fal <- predict.glm(f.a.l, type = "response")

  return(list(fal = fal))
}

#' Cox proportional hazards model for eta-free terms in the (S)AIWKME estimator.
#' @param datalist A list used to calculate the (S)AIWKME estimator including treatment named \code{a}, observed time named \code{obs.t}, censoring indicator (0, censored) named \code{delta}, and baseline covariates used to assign treatment named \code{l}. Notice that all the data in the datalist should be ordered by observed time.
#' @param t0 A predetermined t.
#' @import stats
#' @import survival
#' @return A list including the augmented terms in the numerator with treatment all to 1 named \code{gamma.num.1} and all to 0 named \code{gamma.num.0} and in the denominator with treatment all to 1 named \code{gamma.den.1} and all to 0 named \code{gamma.den.0}; \code{gamma.num.1} and the others are matrix with ordered observed time as rows and patients as columns. More details can be found in references.
#' @export
#' @details More details can be found in references, \code{\link[otrKM]{AIWKME}}, and \code{\link[otrKM]{Genetic.optim.AIWKME}}.
#' @references 
#' {Jiang, R., Lu, W., Song, R., and Davidian, M. (2017) On estimation of optimal treatment regimes for maximizing t‐year survival probability. \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, \bold{79:} 1165-1185. DOI:10.1111/rssb.12201}
Fprep.AIWKME <- function(datalist, t0) {
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
  sur.C <- c(1,rep(fit.C$surv, times = (fit.C$n.event + fit.C$n.censor)))[1:Nevent]

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

#' Logistic regression for observed treatment and instrument used for the (S)IWKMEIV estimator. 
#' @param datalist A list used to calculate the (S)IWKMEIV estimator including treatment named \code{a}, observed time named \code{obs.t}, censoring indicator (0, censored) named \code{delta}, and baseline covariates used to assign treatment named \code{l}. Notice that all the data in the datalist should be ordered by observed time.
#' @param t0 A predetermined t.
#' @import stats
#' @return A list including the probability of receiving instrument given baseline covariates named \code{fzl}, the difference between fal1 and fal0 named \code{deltal}, where fal0 denotes the probability of receiving treatment given baseline covariates and instrument equaling 0, and fal1 denotes the probability of receiving treatment given baseline covariates and instrument equaling 1, and the censoring survival function given baseline covariates and treatment 1 or 0 named \code{surv.C.1} or \code{surv.C.0}.
#' @details More details can be found in references, \code{\link[otrKM]{IWKMEIV}}, and \code{\link[otrKM]{Genetic.optim.IWKMEIV}}.
#' @references 
#' {Xia, J., Zhan, Z., Zhang, J. (2022) Estimating optimal treatment regime in survival contexts using an instrumental variable. Under Review.}
#' @export
Fps.IWKMEIV <- function(datalist, t0) {
  l <- datalist$l
  z <- datalist$z
  a <- datalist$a
  obs.t<-datalist$obs.t
  delta<-datalist$delta

  Nevent=sum(obs.t<=t0)
  rank=order(obs.t)

  # estimate the hazards function
  cox <- coxph(Surv(obs.t, 1-delta)~l+a+z,ties='breslow')

  # estimate baseline hazard function
  cox.a=c(exp(cbind(l,a,z)%*%cox$coefficients))
  cox.a=cox.a[rank]
  delta.rank=(1-delta)[rank]
  den.cox <- sum(cox.a) - cumsum(cox.a) + cox.a
  Lambda0 <- cumsum(delta.rank / den.cox)[1:Nevent] # baseline cumulative hazard function
  # bhest <- c(Lambda0[1], diff(Lambda0))[1:Nevent] # baseline hazard function

  # In the following, the matrix will be ordered observed t times patients
  cox.predict.1<-c(exp(cbind(l,1,z)%*%cox$coefficients))
  cox.predict.0<-c(exp(cbind(l,0,z)%*%cox$coefficients))

  # the first 1 represent a=1, while the second 1 represent z=1
  surv.C.1=exp(-Lambda0%o%cox.predict.1)
  surv.C.0=exp(-Lambda0%o%cox.predict.0)

  f.z.l <- glm(z ~ l, family = binomial(link = "logit"))
  fzl <- predict.glm(f.z.l, type = "response")

  f.a.lz <- glm(a ~ l + z, family = binomial(link = "logit"))
  fal1 <- predict.glm(f.a.lz, newdata = data.frame(l, z = 1), type = "response")
  fal0 <- predict.glm(f.a.lz, newdata = data.frame(l, z = 0), type = "response")
  deltal <- fal1 - fal0

  return(list(surv.C.1=surv.C.1, surv.C.0=surv.C.0, fzl = fzl, deltal = deltal, fal0 = fal0, fal1=fal1))
}

#' Logistic regression for observed treatment and instrument used for the (S)DRKMEIV estimator.
#' @param datalist A list used to calculate the (S)DRKMEIV estimator including treatment named \code{a}, observed time named \code{obs.t}, censoring indicator (0, censored) named \code{delta}, and baseline covariates used to assign treatment named \code{l}. Notice that all the data in the datalist should be ordered by observed time.
#' @param t0 A predetermined t.
#' @import stats
#' @return A list including the probability of receiving instrument given baseline covariates named \code{fzl}, the probability of receiving treatment given baseline covariates and instrument equaling 0 named \code{fal0}, the probability of receiving treatment given baseline covariates and instrument equaling 1 named \code{fal1}, the difference between \code{fal1} and \code{fal0} named \code{deltal}, and the censoring survival function given baseline covariates and treatment 1 or 0 named \code{surv.C.1} or \code{surv.C.0}.
#' @export
#' @details More details can be found in references, \code{\link[otrKM]{DRKMEIV}}, and \code{\link[otrKM]{Genetic.optim.DRKMEIV}}.
#' @references 
#' {Xia, J., Zhan, Z., Zhang, J. (2022) Estimating optimal treatment regime in survival contexts using an instrumental variable. Under Review.}
Fps.DRKMEIV <- function(datalist, t0) {
  l <- datalist$l
  z <- datalist$z
  a <- datalist$a
  obs.t<-datalist$obs.t
  delta<-datalist$delta

  Nevent=sum(obs.t<=t0)
  rank=order(obs.t)

  # estimate the hazards function
  cox.C <- coxph(Surv(obs.t, 1-delta)~l+a+z,ties='breslow')

  # estimate baseline hazard function
  cox.C.a=c(exp(cbind(l,a,z)%*%cox.C$coefficients))
  cox.C.a=cox.C.a[rank]
  delta.C.rank=(1-delta)[rank]
  den.C.cox <- sum(cox.C.a) - cumsum(cox.C.a) + cox.C.a
  Lambda0.C <- cumsum(delta.C.rank / den.C.cox)[1:Nevent] # baseline cumulative hazard function

  # In the following, the matrix will be ordered observed t times patients
  cox.predict.C.1<-c(exp(cbind(l,1,z)%*%cox.C$coefficients))
  cox.predict.C.0<-c(exp(cbind(l,0,z)%*%cox.C$coefficients))

  # the first 1 represent a=1, while the second 1 represent z=1
  surv.C.1=exp(-Lambda0.C%o%cox.predict.C.1)
  surv.C.0=exp(-Lambda0.C%o%cox.predict.C.0)

  f.z.l <- glm(z ~ l, family = binomial(link = "logit"))
  fzl <- predict.glm(f.z.l, type = "response")

  f.a.lz <- glm(a ~ l + z, family = binomial(link = "logit"))
  fal1 <- predict.glm(f.a.lz, newdata = data.frame(l, z = 1), type = "response")
  fal0 <- predict.glm(f.a.lz, newdata = data.frame(l, z = 0), type = "response")
  deltal <- fal1 - fal0

  return(list(surv.C.1=surv.C.1, surv.C.0=surv.C.0, fzl = fzl, deltal = deltal, fal0 = fal0, fal1=fal1))
}

#' Cox proportional hazards model for eta-free terms in the (S)DRKMEIV estimator.
#' @param datalist A list used to calculate the (S)DRKMEIV estimator including treatment named \code{a}, observed time named \code{obs.t}, censoring indicator (0, censored) named \code{delta}, and baseline covariates used to assign treatment named \code{l}. Notice that all the data in the datalist should be ordered by observed time.
#' @param ps A list including the probability of receiving instrument given baseline covariates named \code{fzl}, the probability of receiving treatment given baseline covariates and instrument equaling 0 named \code{fal0}, the probability of receiving treatment given baseline covariates and instrument equaling 1 named \code{fal1}, and the difference between fal1 and fal0 named \code{deltal}. \code{\link[otrKM]{Fps.DRKMEIV}} can produce \code{ps} by positing logistic model.
#' @param t0 A predetermined t.
#' 
#' @import stats
#' @import survival
#' @return A list including estimates \eqn{\hat{\gamma}_1(\boldsymbol{L};s)}{hat.gamma.1(L;s)} with treatment all to 1 named \code{gamma.num.1} and all to 0 named \code{gamma.num.0}, \eqn{\hat{\gamma}_1'(\boldsymbol{L};s)}{hat.gamma.1'(L;s)} with treatment all to 1 named \code{gammaa.num.1} and all to 0 named \code{gammaa.num.0}, \eqn{\hat{\gamma}_2(\boldsymbol{L};s)}{hat.gamma.2(L;s)} with treatment all to 1 named \code{gamma.den.1} and all to 0 named \code{gamma.den.0}, and \eqn{\hat{\gamma}_2'(\boldsymbol{L};s)}{hat.gamma.2'(L;s)} with treatment all to 1 named \code{gammaa.den.1} and all to 0 named \code{gammaa.den.0}; \code{gamma.num.1} and the others are matrix with ordered observed time as rows and patients as columns. There are also estimates for the last term of the (S)DRIWKMEIV estimator. More details can be found in references.
#' 
#' @export
#' @details More details can be found in references, \code{\link[otrKM]{DRKMEIV}}, and \code{\link[otrKM]{Genetic.optim.DRKMEIV}}.
#' @references 
#' {Xia, J., Zhan, Z., Zhang, J. (2022) Estimating optimal treatment regime in survival contexts using an instrumental variable. Under Review.}
Fprep.DRKMEIV <- function(datalist, ps, t0){
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
  n=length(a)

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

  gammaa.num.1=t(t(surv.T.10*(bhest%o%cox.predict.10))*fal0)

  gamma.num.1=t(t(surv.T.11*(bhest%o%cox.predict.11))*fal1/sapply(deltal,clipp))-t(t(surv.T.10*(bhest%o%cox.predict.10))*fal0/sapply(deltal,clipp))

  gammaa.den.1=t(t(surv.T.10)*fal0)

  gamma.den.1=t(t(surv.T.11)*fal1/sapply(deltal,clipp))-t(t(surv.T.10)*fal0/sapply(deltal,clipp))

  cox.predict.01<-c(exp(cbind(l,0,1,l*0)%*%cox$coefficients))
  cox.predict.00<-c(exp(cbind(l,0,0,l*0)%*%cox$coefficients))

  surv.T.01=exp(-Lambda0%o%cox.predict.01)
  surv.T.00=exp(-Lambda0%o%cox.predict.00)

  gammaa.num.0=-t(t(surv.T.00*(bhest%o%cox.predict.00))*(1-fal0))

  gamma.num.0=-t(t(surv.T.01*(bhest%o%cox.predict.01))*(1-fal1)/sapply(deltal,clipp))+t(t(surv.T.00*(bhest%o%cox.predict.00))*(1-fal0)/sapply(deltal,clipp))

  gammaa.den.0=-t(t(surv.T.00)*(1-fal0))

  gamma.den.0=-t(t(surv.T.01)*(1-fal1)/sapply(deltal,clipp))+t(t(surv.T.00)*(1-fal0)/sapply(deltal,clipp))
  
  # estimate the censoring hazards function
  cox.C <- coxph(Surv(obs.t, 1-delta)~l+a+z,ties='breslow')

  # estimate baseline hazard function
  cox.C.a=c(exp(cbind(l,a,z)%*%cox.C$coefficients))[rank]
  delta.C.rank=(1-delta)[rank]
  den.C.cox <- sum(cox.C.a) - cumsum(cox.C.a) + cox.C.a
  Lambda0.C <- cumsum(delta.C.rank / den.C.cox)[1:Nevent] # baseline cumulative hazard function
  bhest.C <- c(Lambda0.C[1], diff(Lambda0.C))[1:Nevent] # baseline hazard function

  # In the following, the matrix will be ordered observed t times patients
  cox.predict.C.1<-c(exp(cbind(l,1,z)%*%cox.C$coefficients))
  cox.predict.C.0<-c(exp(cbind(l,0,z)%*%cox.C$coefficients))

  # the first 1 represent a=1, while the second 1 represent z=1
  surv.C.1=exp(-Lambda0.C%o%cox.predict.C.1)
  surv.C.0=exp(-Lambda0.C%o%cox.predict.C.0)

  dNC=cbind(diag((delta.C.rank)[1:Nevent]), matrix(0, nrow = Nevent, ncol =n-Nevent))
  Y=dNC
  Y[upper.tri(Y,diag=TRUE)] = 1
  dMC1= dNC- bhest.C%o%cox.predict.C.1*Y
  dMC0 =dNC- bhest.C%o%cox.predict.C.0*Y

  surv.T.1=surv.T.11*z+surv.T.10*(1-z)
  surv.T.0=surv.T.01*z+surv.T.00*(1-z)

  temp1=rbind(0,apply(dMC1/rbind(surv.C.1[-1,],1)/surv.T.1,2,cumsum)[-Nevent,])
  temp0=rbind(0,apply(dMC0/rbind(surv.C.0[-1,],1)/surv.T.0,2,cumsum)[-Nevent,])

  mat.num.1=(bhest%o%cox.predict.11*z+bhest%o%cox.predict.10*(1-z))*surv.T.1*temp1
  mat.num.0=(bhest%o%cox.predict.01*z+bhest%o%cox.predict.00*(1-z))*surv.T.0*temp0

  mat.den.1=surv.T.1*temp1
  mat.den.0=surv.T.0*temp0

  return(list(  gammaa.num.1=gammaa.num.1,
                gammaa.num.0=gammaa.num.0,
                gamma.num.1=gamma.num.1,
                gamma.num.0=gamma.num.0,
                gammaa.den.1=gammaa.den.1,
                gammaa.den.0=gammaa.den.0,
                gamma.den.1=gamma.den.1,
                gamma.den.0=gamma.den.0,
                mat.num.1=mat.num.1,
                mat.num.0=mat.num.0,
                mat.den.1=mat.den.1,
                mat.den.0=mat.den.0))
}