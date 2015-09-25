### parameter solution by numerical approach. Failed (inaccurate & too expensive)

library(cubature)   #for hypercubed integration

### mv.normal pdf
dmvnorm <- function (x, mean, sigma, log = FALSE) {
  if (is.vector(x)) {
    x <- matrix(x, ncol = length(x))
  }
  if (missing(mean)) {
    mean <- rep(0, length = ncol(x))
  }
  if (missing(sigma)) {
    sigma <- diag(ncol(x))
  }
  if (NCOL(x) != NCOL(sigma)) {
    stop("x and sigma have non-conforming size")
  }
  if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps),
                   check.attributes = FALSE)) {
    stop("sigma must be a symmetric matrix")
  }
  if (length(mean) != NROW(sigma)) {
    stop("mean and sigma have non-conforming size")
  }
  distval <- mahalanobis(x, center = mean, cov = sigma)
  logdet <- sum(log(eigen(sigma, symmetric = TRUE, only.values = TRUE)$values))
  logretval <- -(ncol(x) * log(2 * pi) + logdet + distval)/2
  if (log)
    return(logretval)
  exp(logretval)
}
adaptIntegrate(f=dmvnorm, mean=rep(0,5), sigma=Sigma, lowerLimit=rep(-100,5), upperLimit=rep(100,5), 
               doChecking=F, maxEval=1000000)

### E(Vi|D=1)
vPhiphi = function (x, alpha0, alpha1, mean, sigma) {
  Phi = (1+ exp(alpha0 + x%*%alpha1)^-1)^-1
  Phi0 = (1+ exp(alpha0)^-1)^-1
  phi = dmvnorm(x,mean, sigma)
  return(x*Phi*phi/Phi0)
}

vPhiphi0 = function (x, alpha0, alpha1, mean, sigma) {
  Phi = 1 - (1+ exp(alpha0 + x%*%alpha1)^-1)^-1
  Phi0 = 1 - (1+ exp(alpha0)^-1)^-1
  phi = dmvnorm(x,mean, sigma)
  return(x*Phi*phi/Phi0)
}

adaptIntegrate(vPhiphi0, fDim=5, lowerLimit=rep(-100,5), upperLimit=rep(100,5), alpha0=0, alpha1=rep(1,5), mean=rep(0,5), sigma=Sigma, 
                      maxEval=2000000)
