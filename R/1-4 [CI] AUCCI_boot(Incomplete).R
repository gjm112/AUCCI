### Part I. CI functions
## 1-4 [CI] AUCCI_boot(Incomplete):  1.4.1 AUCCI.boot(verification bias)

## 1.4.1 AUCCI.boot ##################################################################### 
### Bootstrap CI for AUC with verification bias
AUCCI.boot = function(x,R,alpha,fun,logit=TRUE,...) {
  # x: dataframe of the sample data
  # R: number of replicates
  # fun: AUC function
  n = dim(x)[1]
  boot.sample.indices = matrix(sample(1:n,size=n*R,replace=TRUE), R, n)
  
  AUC.hats = rep(NA,R)
  for (i in 1:R) {AUC.hats[i] = fun(x[boot.sample.indices[i,],],...)}
  mu.hat = fun(x,...)
  result = list(boot.statistics = AUC.hats)
  if (logit == TRUE) {
    mu.hat.logit = logit(mu.hat)
    AUC.hats.logit = logit(AUC.hats)
    se.logit = sd(AUC.hats.logit, na.rm=TRUE)
    interval.logit = mu.hat.logit + c(-1,1)*(-qnorm(alpha/2))*se.logit
    interval = (exp(interval.logit)^-1 +1)^-1
    result$interval = interval
    result$se.logit = se.logit
  } else {
    se = sd(AUC.hats, na.rm=TRUE)
    interval = mu.hat + c(-1,1)*(-qnorm(alpha/2))*se
    result$interval = interval
    result$se = se
  }
  print( interval )
  return( result )
}

## Example ############################################################################## 
AUCCI.boot(temp,R=100,alpha=.05,fun=AUC.verif,method="BG")