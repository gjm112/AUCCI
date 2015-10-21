### Part I. CI functions
## 1-4 [CI] AUCCI_boot(Incomplete):  1.4.1 AUCCI.boot(verification bias)

## 1.4.1 AUCCI.boot ##################################################################### 
### Bootstrap CI for AUC with verification bias
AUCCI.boot = function(x,R,alpha,fun,print.boots = FALSE, variance=FALSE, LT=TRUE,...) {
  # x: dataframe of the sample data
  # R: number of replicates
  # fun: AUC function
  n = dim(x)[1]
  boot.sample.indices = matrix(sample(1:n,size=n*R,replace=TRUE), R, n)
  AUC.hat = fun(x,...)
  AUC.hats = rep(NA,R)
  for (i in 1:R) {AUC.hats[i] = fun(x[boot.sample.indices[i,],],...)}
  mu.hat = fun(x,...)
  result = list(AUC.hat = AUC.hat)
  if (print.boots) {result$boot.statistics = AUC.hats}
  if (LT == TRUE) {
    mu.hat.LT = logit(mu.hat)
    AUC.hats.LT = logit(AUC.hats)
    se.LT = sd(AUC.hats.LT, na.rm=TRUE)
    interval.LT = mu.hat.LT + c(-1,1)*(-qnorm(alpha/2))*se.LT
    interval = (exp(interval.LT)^-1 +1)^-1
    result$CI = interval
    if (variance) {result$variance = (se.LT*AUC.hat*(1-AUC.hat))^2}
  } else {
    se = sd(AUC.hats, na.rm=TRUE)
    interval = mu.hat + c(-1,1)*(-qnorm(alpha/2))*se
    result$CI = interval
    if (variance) {result$variance = se^2}
  }
  return( result )
}

## Example ############################################################################## 
AUCCI.boot(temp,R=100,alpha=.05,fun=AUC.verif,method="BG")
