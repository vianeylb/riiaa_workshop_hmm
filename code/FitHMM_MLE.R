## Code from Zuchinni's HMMs for time series book

#Transform sigma by taking log
norm_HMM_pn2pw <- function(m,mu,sigma,gamma,delta=NULL,stationary=TRUE)
{
  tsigma <- log(sigma)
  foo <- log(gamma/diag(gamma))
  tgamma <- as.vector(foo[!diag(m)])
  if(stationary) 
  {tdelta <-NULL} 
  else 
  {tdelta<-log(delta[-1]/delta[1])}
  parvect <- c(mu,tsigma,tgamma,tdelta)
  return(parvect)
}

# Transform normal working parameters to natural parameters
norm_HMM_pw2pn <- function(m,parvect,stationary=TRUE)
{
  #mu is untransformed
  mu <- parvect[1:m]
  #Take exponent of tsigma
  sigma <- exp(parvect[(m+1):(2*m)])
  gamma <- diag(m)
  gamma[!gamma] <- exp(parvect[(2*m+1):(m+m*m)])
  gamma <- gamma/apply(gamma,1,sum)
  if(stationary) 
  {delta<-solve(t(diag(m)-gamma+1),rep(1,m))} 
  else
  {foo<-c(1,exp(parvect[(m+m*m+1):(m*m+2*m-1)]))
  delta<-foo/sum(foo)}
  return(list(mu=mu,sigma=sigma,gamma=gamma,delta=delta))
}

#Computing minus the log-likelihood from the working parameters
norm_HMM_mllk <- function(parvect, x, m, stationary=TRUE)
{
  if(m==1) 
    return (-sum(dnorm(x, mean=parvect[1], sd=exp(parvect[2]), log=TRUE)))
  n <- length(x)
  pn <- norm_HMM_pw2pn(m, parvect, stationary=stationary)
  foo <- pn$delta*dnorm(x[1], pn$mu, pn$sigma)
  sumfoo <- sum(foo)
  lscale <- log(sumfoo)
  foo <- foo/sumfoo
  for (i in 2:n)
  {
    if (!is.na(x[i]))
    {P<-dnorm(x[i], pn$mu, pn$sigma)}
    else
    {P<- rep(1,m)}
    foo <- foo%*%pn$gamma*P
    sumfoo <- sum(foo)
    lscale <- lscale + log(sumfoo)
    foo <- foo/sumfoo
  }
  mllk <- -lscale
  return(mllk)
}

#Computing MLE from natural parameters
norm_HMM_mle <- function(x, m, mu0, sigma0, gamma0, delta0=NULL, stationary=TRUE)
{
  parvect0 <- norm_HMM_pn2pw(m, mu0, sigma0, gamma0, delta0, stationary=stationary)
  mod <- nlm(norm_HMM_mllk, parvect0, x=x,m=m, stationary=stationary, print.level =2)
  pn <- norm_HMM_pw2pn(m=m, mod$estimate, stationary=stationary)
  mllk <- mod$minimum

  return(list(pn, mllk))
}

## Run HMM ----------------------------------

m <- 2
mu0 <- c(4.5, 9)
sigma0 <- c(3, 3)
gamma0 <- matrix(c(0.8, 0.2, 
                  0.2, 0.8), nrow = m, byrow = TRUE)
delta0 <- c(0.5, 0.5)
mod <- list(m=m, mu=mu0, sigma=sigma0, gamma=gamma0, delta=delta0)
norm_mle <- norm_HMM_mle(obs, m, mu0, sigma0, gamma0, delta0, stationary=TRUE)
