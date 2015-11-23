### Part I. CI functions
## 1-4 [CI] AUCCI_boot(Incomplete):  1.4.1 AUCCI.boot(verification bias)

## 1.4.1 AUCCI.boot ##################################################################### 
### Bootstrap CI for AUC with verification bias
# norm.inter from the package boot
norm.inter <- function (t, alpha) {
  if (is.list(t)) {t <- t[sapply(t, is.finite)]} else{t <- t[is.finite(t)]}
  R <- length(t)
  rk <- (R + 1) * alpha
  if (is.na(rk)) {return(matrix(NA,2,2))}   # debugging
  if (!all(rk > 1 & rk < R)) 
    warning("extreme order statistics used as endpoints")
  k <- trunc(rk)
  inds <- seq_along(k)
  out <- inds
  kvs <- k[k > 0 & k < R]
  tstar <- sort(t, partial = sort(union(c(1, R), c(kvs, kvs + 
                                                     1))))
  ints <- (k == rk)
  if (any(ints)) 
    out[inds[ints]] <- tstar[k[inds[ints]]]
  out[k == 0] <- tstar[1L]
  out[k == R] <- tstar[R]
  not <- function(v) xor(rep(TRUE, length(v)), v)
  temp <- inds[not(ints) & k != 0 & k != R]
  temp1 <- qnorm(alpha[temp])
  temp2 <- qnorm(k[temp]/(R + 1))
  temp3 <- qnorm((k[temp] + 1)/(R + 1))
  tk <- tstar[k[temp]]
  tk1 <- tstar[k[temp] + 1L]
  out[temp] <- tk + (temp1 - temp2)/(temp3 - temp2) * (tk1 - tk)
  cbind(round(rk, 2), out)
}

AUCCI.boot = function(data,R,alpha=0.05,base.fun=AUC.verif,CI.method, print.boots = FALSE, variance=FALSE, LT=TRUE, type="portrait", ...) {
  # x: dataframe of the sample data
  # R: number of replicates
  # base.fun: AUC function
  len = length(CI.method)
  n = dim(data)[1]
  boot.sample.indices = matrix(sample(1:n,size=n*R,replace=TRUE), R, n)
  AUC.hats <- data.frame(matrix(NA,len,R))  # Bootstrapped AUC hats
  CIs <- data.frame(matrix(NA,len,5))     #5 = AUC.hat, Wald lb, ub, BCA lb, ub
  colnames(CIs) = c("AUC.hat", "Wald.lb", "Wald.ub", "BCA.lb", "BCA.ub")
  CIs[,1] <- AUC.hat <- t(base.fun(data, CI.method=CI.method, ...))
  for (j in 1:R) {AUC.hats[,j] = t(base.fun(data[boot.sample.indices[j,],], CI.method=CI.method, ...))}
  AUC.hat.LT = logit(CIs[,1], LT = LT)
  AUC.hats.LT = logit(AUC.hats, LT = LT)
  se.LT = apply(AUC.hats.LT,1,sd, na.rm=TRUE)
  
  ### 1. Wald CI's
  Wlb.LT = AUC.hat.LT -(-qnorm(alpha/2))*se.LT; Wlb = expit(Wlb.LT, LT = LT)
  Wub.LT = AUC.hat.LT +(-qnorm(alpha/2))*se.LT; Wub = expit(Wub.LT, LT = LT)
  
  CIs[,2] = Wlb
  CIs[,3] = Wub  
  
  ### 2. BCA CI's - codes excerpted from the package boot
  # For BCA, LT is ignored
  # rate of (t < t0)  where t: B/S estimates , t0 : original est'
  t.r <- apply(AUC.hats < matrix(rep(CIs[,1], R), len, R), 1, mean, na.rm = TRUE)
  w <- qnorm(t.r)
  zalpha = qnorm(c(0,1) + alpha * c(1,-1) /2)
    
  # Jackknife stats
  AUC.hats.J <- data.frame(matrix(NA,len,n))
  for (j in 1:n) {
    AUC.hats.J[,j] <- t(base.fun(data[-j,], CI.method=CI.method, ...))
  }
  mean.J = apply(AUC.hats.J, 1, mean, na.rm=TRUE)
  L = AUC.hats.J - matrix(rep(mean.J, n), len, n)
  a <- apply(L^3,1,sum) / apply(L^2, 1, sum)^1.5
  a.adj.1 <- pnorm(w + (w + zalpha[1])/(1 - a * (w + zalpha[1])))
  a.adj.2 <- pnorm(w + (w + zalpha[2])/(1 - a * (w + zalpha[2])))
  qq <- list()
  for (i in 1:len) {
    qq[[i]] <- norm.inter(as.matrix(AUC.hats[i,]), c(a.adj.1[i], a.adj.2[i]))
  }
  for (i in 1:len)  {CIs[i,4] <- qq[[i]][1,2L]; CIs[i,5] <- qq[[i]][2,2L]}
  rownames(CIs) <- CI.method
  
  CIs = CIs[,c(1,4,5,2,3)]
  if (type=="portrait") {
    result = list(CI = CIs)
  } else if (type=="landscape1"|type=="landscape2") {
    CI.BCA = as.data.frame(t(CIs[,c(2,3)]))
    colnames(CI.BCA) = CI.method
    rownames(CI.BCA) = c("BCA.lb", "BCA.ub")
    
    CI.Wald = as.data.frame(t(CIs[,c(4,5)]))
    colnames(CI.Wald) = CI.method
    rownames(CI.Wald) = c("Wald.lb", "Wald.ub")
    
    if (type=="landscape2") {
      CI.BCA = t(as.vector(as.matrix(CI.BCA)))
      CI.Wald = t(as.vector(as.matrix(CI.Wald)))
      colnames(CI.BCA) = c(paste0(rep(CI.method,each=2),c(".lb",".ub")))
      colnames(CI.Wald) = c(paste0(rep(CI.method,each=2),c(".lb",".ub")))
      rownames(CI.BCA) = "BCA"
      rownames(CI.Wald) = "Wald"}
    
    result = list(CI = CI.BCA, CI.Wald = CI.Wald, AUC.hat = AUC.hat)
    
  } else {
    result = list(CI = CIs)
    warning ("Wrong type of print: portrait type is used instead.")
  }
  
  if (variance) {result$variance = (se.LT*ifelse(LT, AUC.hat*(1-AUC.hat), 1))^2}
  if (print.boots) {result$boot.statistics = AUC.hats}
  return( result )
}