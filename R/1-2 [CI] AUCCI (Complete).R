### Part I. CI functions
## 1-2 [CI] AUCCI(Complete): 1.2.1 AUCCI, 1.2.2 CI.i

## 1.2.1 AUC CI #######################################################################
AUCCI <- function(data, method, disease="disease", marker="marker", alpha=0.05, LT=FALSE, variance = FALSE, cc = FALSE) {
  ### note for developers
  ## method
  # "HanleyMcNeilWald", "HanleyMcNeilExponential", "HanleyMcNeilScore", ...
  ## variance: logical
  # to show estimated variance also or not?
  ## cc: logical
  # continuity correction (for Wilson Score only)
  ## LT: logit transformation (HMW,HME, NE, HMS, NS, CM, Bb, MW, DL)  ( RG, HalpM, DB, DG, CP, WS)
  # For Score-like methods(WS, HalpM) or parametric methods(RG,DB,DG), LT is not necessary.
  # returning variance when logit=T is var(logit(theta))
  
  # basic stats
  x = data[data[,disease]==0, marker]
  y = data[data[,disease]==1, marker]
  n.x = length(x) ; n.y = length(y) ; n = n.x + n.y
  A = AUC(x, y, n.x, n.y)
  A.LT = logit(A, LT=LT)
  
  if (method == "HanleyMcNeilWald") {
    Q = Q.stat(x, y, n.x, n.y); Q1 = Q$Q1; Q2 = Q$Q2
    V.LT = V.HM(A, n.x, n.y, Q1, Q2,LT=LT)
    CI.LT = CI.base(A.LT,V.LT,alpha)
  } 
  
  else if (method == "HanleyMcNeilExponential") {
    V.LT = V.HM(A, n.x, n.y, LT=LT)
    CI.LT = CI.base(A.LT,V.LT,alpha)
  }
  
  else if (method == "HanleyMcNeilScore") {
    CI = multiroot(HMS.equation, c(0.5,0.9), AUC.hat=A, n.x=n.x, n.y=n.y, alpha=alpha, LT=LT)$root  #original scale
    CI.LT = logit(CI, LT=LT)
    if (variance) {V.LT = V.HM(A, n.x, n.y, LT=LT)}
  }
  
  else if (method == "NewcombeExponential") {
    V.LT = V.NC(A, n.x, n.y, LT=LT)
    if (V.LT < 0) { warning("Negative estimated variance set to zero")
                 V.LT = 0 }
    CI.LT = CI.base(A.LT,V.LT,alpha)
  }
  
  else if (method == "NewcombeScore") {
    CI = multiroot(NC.equation, c(0.5,0.9), AUC.hat=A, n.x=n.x, n.y=n.y, alpha=alpha, LT=LT)$root
    CI.LT = logit(CI, LT=LT)
    if (variance) {
      V.LT = V.NC(A, n.x, n.y, LT=LT)
      if (V.LT < 0) { warning("Negative estimated variance set to zero")
                   V.LT = 0 }      
    }
  }
  
  else if (method == "CortesMohri") {
    V.LT = V.CM(A, n.x, n.y, LT=LT)
    CI.LT = CI.base(A.LT,V.LT,alpha)
  }
  
  else if (method == "ReiserGuttman") {
    probit = probit.RG(x,y,n.x,n.y)
    d = probit$delta; V = probit$V
    A = pnorm(d)
    d.CI = d + c(-1,1)*qnorm(1-alpha/2)*sqrt(V)
    CI = pnorm(d.CI)
    CI.LT = logit(CI, LT=LT) # logit transformation not available for RG (only for presentation purpose)
    V.LT = NA                # variance not available for RG(The V above is variance in probit scale)
  }
  
  else if (method == "Bamber") {
    V.LT = V.Bm(x, y, n.x, n.y, A, LT=LT)
    CI.LT = CI.base(A.LT,V.LT,alpha)
  }

  else if (method == "HalperinMee") {
    CI = multiroot(Halperin.equation, c(0.5,0.9), AUC.hat=A, x=x, y=y, n.x=n.x, n.y=n.y, alpha=alpha, LT=LT)$root
    CI.LT = logit(CI, LT=LT)
    if (variance) {
      V.LT=V.Halperin(x, y, n.x=n.x, n.y=n.y, LT=LT) 
    }
  }
  
  else if (method == "DoubleBeta") {     #variance(or alp) has an unstable solution
    CI = DB(A, n.x, n.y, alpha, LT=LT)
    CI.LT = logit(CI, LT=LT)
    if (variance) {
      alp = multiroot(DB.equation, c(1,3), AUC.hat=A, n.x=n.x, n.y=n.y, alpha=1,LT=LT)$root[1]
      V.LT = V.DB(alp, n.x, n.y, LT=LT)
    }
  }
  
  else if (method == "DoubleGaussian") {  #variance(or alp) has an unstable solution
    CI = DG(A, n.x, n.y, alpha, LT=LT)
    CI.LT = logit(CI, LT=LT)
    if (variance) {
      delta = multiroot(DG.equation, c(1,3), AUC.hat=A, n.x=n.x, n.y=n.y, alpha=1,LT=LT)$root[1]
      V.LT = V.DG(delta, n.x, n.y, LT=LT)
    }
  }
  
  else if (method == "MannWhitney") {
    V.LT = S.stat(x, y, n.x, n.y, A, LT=LT)^2 * n / (n.x * n.y)
    CI.LT = CI.base(A.LT,V.LT,alpha)
  }
  
  else if (method == "DeLong") {
    V.LT = V.DeLong(x, y, n.x, n.y, A, LT=LT)
    CI.LT = CI.base(A.LT,V.LT,alpha)
  }  
  
  ##### WilsonScore and ClopperPearson methods are excluded from this study ##########
  else if (method == "WilsonScore") {
    z.a2 = qnorm(1-alpha/2)
    CI = (2*n*A + z.a2^2 +cc*c(-1,1) +c(-1,1)* z.a2* sqrt(z.a2^2 + cc*c(-2,2) - 1/n*cc + 4*A*(n*(1-A) +1*cc)))/(2*(n+z.a2^2))
  }
  
  # V for CP???
  else if (method == "ClopperPearson") {
    n = n.x + n.y; k = round(A * n)
    f1 = qf(alpha/2,2*k,2*(n-k+1))
    f2 = qf(1-alpha/2,2*(k+1),2*(n-k))
    ##  V = ??
    CI = c( (k*f1)/(n-k+1+k*f1), ((k+1)*f2)/(n-k+(k+1)*f2) )
  }
  ##### WilsonScore and ClopperPearson methods are excluded from this study ##########
  
  result = list(AUC.hat = A, CI = expit(CI.LT,LT=LT));  if (variance) {result$V.hat = V.LT}
  return(result)
}


## 1.2.2 lumpsum AUCCI for simulation purpose #########################################
# CI's for individual elements of a sample (1. landscape form)  2 (lb,ub) x 17
CI.i <- function(data, fun=AUCCI, CI.method=CI.methods, type="landscape2", ...) {
  if (type == "landscape1") {
    CI.i = as.data.frame(matrix(NA,2,length(CI.method)+1))
    names(CI.i) = c("AUC.hat", CI.method)
    CI.i$AUC.hat = fun(data, CI.method[1],...)$AUC.hat
    for (i in 1:length(CI.method)) { CI.i[,i+1] = fun(data, CI.method = CI.method[i], ...)$CI }
  } else if (type == "landscape2") {
    CI.i = as.data.frame(matrix(NA,1,length(CI.method)*2+1))
    names(CI.i) = c("AUC.hat",paste0(rep(CI.methods,each=2),c(".lb",".ub")))
    CI.i$AUC.hat = fun(data, CI.method[1],...)$AUC.hat
    for (i in 1:length(CI.method)) { CI.i[1,(2*i):(2*i+1)] = fun(data, CI.method = CI.method[i], ...)$CI }  
  } else if (type=="portrait") {
    CI.i = as.data.frame(matrix(NA,length(CI.method)+1,2))
    names(CI.i) = c("lowerbound", "upperbound")
    rownames(CI.i)  = c("AUC.hat", CI.method)
    CI.i[1,] = fun(data, CI.method[1],...)$AUC.hat
    for (i in 1:length(CI.method)) { CI.i[i+1,] = fun(data, CI.method = CI.method[i], ...)$CI }  
  } else {stop("Form is none of landscape1, landscape2, or portrait")}    
  return(CI.i)
}
