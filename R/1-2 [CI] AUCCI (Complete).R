### Part I. CI functions
## 1-2 [CI] AUCCI(Complete): 1.2.1 AUCCI, 1.2.2 CI.i

## 1.2.1 AUC CI #######################################################################
AUCCI <- function(data, method, disease="disease", marker="marker", alpha=0.05, variance = FALSE, cc = FALSE) {
  ### note for developers
  ## method
  # "HanleyMcNeilWald", "HanleyMcNeilExponential", "HanleyMcNeilScore", ...
  ## variance: logical
  # to show estimated variance also or not?
  ## cc: logical
  # continuity correction (for Wilson Score only)
  
  # basic stats
  x = data[data[,disease]==0, marker]
  y = data[data[,disease]==1, marker]
  n.x = length(x) ; n.y = length(y) ; n = n.x + n.y
  A = AUC(x, y, n.x, n.y)
  
  if (method == "HanleyMcNeilWald") {
    Q = Q.stat(x, y, n.x, n.y); Q1 = Q$Q1; Q2 = Q$Q2
    V = V.HM(A, n.x, n.y, Q1, Q2)
    CI = CI.base(A,V,alpha)
  } 
  
  else if (method == "HanleyMcNeilExponential") {
    V = V.HM(A, n.x, n.y)
    CI = CI.base(A,V,alpha)  
  }
  
  else if (method == "HanleyMcNeilScore") {
    CI = multiroot(HMS.equation, c(0.01,0.99), AUC.hat=A, n.x=n.x, n.y=n.y, alpha=alpha)$root
    if (variance) {V = V.HM(A, n.x, n.y)}
  }
  
  else if (method == "NewcombeExponential") {
    V = V.NC(A, n.x, n.y)
    if (V < 0) { warning("Negative estimated variance set to zero")
                 V = 0 }
    CI = CI.base(A,V,alpha)  
  }
  
  else if (method == "NewcombeScore") {
    CI = multiroot(NC.equation, c(0.01,0.99), AUC.hat=A, n.x=n.x, n.y=n.y, alpha=alpha)$root
    if (variance) {
      V = V.NC(A, n.x, n.y)
      if (V < 0) { warning("Negative estimated variance set to zero")
                   V = 0 }      
    }
  }
  
  else if (method == "CortesMohri") {
    V = V.CM(A, n.x, n.y)
    CI = CI.base(A,V,alpha)
  }
  
  else if (method == "ReiserGuttman") {
    mu.x = mean(x); mu.y = mean(y)
    sig.x = sd(x);  sig.y = sd(y); s2 = sig.x^2 + sig.y^2
    d = (mu.y-mu.x)/sqrt(s2)
    f = s2^2/((sig.x^4/(n.x-1))+(sig.y^4/(n.y-1)))
    M = s2/((sig.x^2/n.x)+(sig.y^2/n.y))
    V = 1/M + d^2/(2*f)
    d.CI = d + c(-1,1)*qnorm(1-alpha/2)*sqrt(V)
    CI = pnorm(d.CI)
  }
  
  # V for CP???
  else if (method == "ClopperPearson") {
    n = n.x + n.y; k = round(A * n)
    f1 = qf(alpha/2,2*k,2*(n-k+1))
    f2 = qf(1-alpha/2,2*(k+1),2*(n-k))
    ##  V = ??
    CI = c( (k*f1)/(n-k+1+k*f1), ((k+1)*f2)/(n-k+(k+1)*f2) )
  }
  
  else if (method == "Bamber") {
    V = V.Bm(x, y, n.x, n.y, A)
    CI = CI.base(A,V,alpha)
  }
  
  else if (method == "WilsonScore") {
    z.a2 = qnorm(1-alpha/2)
    CI = (2*n*A + z.a2^2 +cc*c(-1,1) +c(-1,1)* z.a2* sqrt(z.a2^2 + cc*c(-2,2) - 1/n*cc + 4*A*(n*(1-A) +1*cc)))/(2*(n+z.a2^2))
  }
  
  else if (method == "HalperinMee") {
    CI = Halperin(x, y, n.x, n.y, alpha)
    if (variance) {
      Halperin.stat = Halperin.stat(x, y, n.x, n.y)
      AUC.hat = Halperin.stat$AUC.hat; N.J.hat = Halperin.stat$N.J.hat
      V = AUC.hat*(1-AUC.hat)/N.J.hat
    }
  }
  
  # V for DB: TBD
  else if (method == "DoubleBeta") {
    CI = DB(A, n.x, n.y, alpha)
  }
  
  # V for DG: TBD
  else if (method == "DoubleGaussian") {
    CI = DG(A, n.x, n.y, alpha)
  }
  
  else if (method == "MannWhitney") {
    V = S.stat(x, y, n.x, n.y)^2 * n / (n.x * n.y)
    CI = CI.base(A,V,alpha)
  }
  
  else if (method == "MannWhitneyLT") {
    V = S.stat(x, y, n.x, n.y)^2 * n / (n.x * n.y)
    CI = CI.base(log(A/(1-A)), V/ (A*(1-A))^2, alpha)
    CI = exp(CI) / (exp(CI) + 1)
  }
  
  else if (method == "DeLong") {
    V = V.DeLong(x, y, n.x, n.y, A)
    CI = CI.base(A,V,alpha)
  }  
  
  result = list(AUC.hat = A, CI = CI);  if (variance) {result$V.hat = V}
  return(result)
}


## 1.2.2 lumpsum AUCCI for simulation purpose #########################################
# CI's for individual elements of a sample (1. landscape form)  2 (lb,ub) x 17
CI.i <- function(data, method=CI.methods, type="landscape2", ...) {
  if (type == "landscape1") {
    CI.i = as.data.frame(matrix(NA,2,length(method)+1))
    names(CI.i) = c("AUC.hat", method)
    CI.i$AUC.hat = AUCCI(data, method[1])$AUC.hat
    for (i in 1:length(method)) { CI.i[,i+1] = AUCCI(data, method = method[i], ...)$CI }
  } else if (type == "landscape2") {
    CI.i = as.data.frame(matrix(NA,1,length(method)*2+1))
    names(CI.i) = c("AUC.hat",paste0(rep(CI.methods,each=2),c(".lb",".ub")))
    CI.i$AUC.hat = AUCCI(data, method[1])$AUC.hat
    for (i in 1:length(method)) { CI.i[1,(2*i):(2*i+1)] = AUCCI(data, method = method[i], ...)$CI }  
  } else if (type=="portrait") {
    CI.i = as.data.frame(matrix(NA,length(method)+1,2))
    names(CI.i) = c("lowerbound", "upperbound")
    rownames(CI.i)  = c("AUC.hat", method)
    CI.i[1,] = AUCCI(data, method[1])$AUC.hat
    for (i in 1:length(method)) { CI.i[i+1,] = AUCCI(data, method = method[i], ...)$CI }  
  } else {stop("Form is none of landscape1, landscape2, or portrait")}    
  return(CI.i)
}
