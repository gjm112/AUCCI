### Part I. CI functions
## 1-2 [CI] AUCCI(Complete): 1.2.1 AUCCI, 1.2.2 CI.i

## 1.2.1 AUC CI #######################################################################
AUCCI <- function(data, CI.method, disease="disease", marker="marker", alpha=0.05, LT=FALSE, variance = FALSE, cc = FALSE, ...) {
  ### note for developers
  ## CI.method: "HanleyMcNeilWald", "HanleyMcNeilExponential", "HanleyMcNeilScore", ...
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
  x = x[!is.na(x)]
  y = y[!is.na(y)]
  n.x = length(x) ; n.y = length(y) ; n = n.x + n.y
  if (n.x == 0) {warning("there is no x(nondiseased)")}
  else if (n.y == 0) {warning("there is no y(diseased)")}
  
  A = AUC(x, y, n.x, n.y, ...)
  A.LT = logit(A, LT=LT)
  start = c( A/2, (A+1)/2)
  start.LT = logit(start, LT=LT)
  
  if (CI.method == "HM1") {
    Q = Q.stat(x, y, n.x, n.y); Q1 = Q$Q1; Q2 = Q$Q2
    pXY = p.XY(x, y, n.x, n.y)
    V.LT = V.HM(A, n.x, n.y, Q1, Q2, pXY, LT=LT,...)
    CI.LT = CI.base(A.LT,V.LT,alpha)
  } 
  
  else if (CI.method == "HM2") {
    pXY = p.XY(x, y, n.x, n.y)
    V.LT = V.HM(A, n.x, n.y, pXY=pXY, LT=LT,...)
    CI.LT = CI.base(A.LT,V.LT,alpha)
  }
  
  else if (CI.method == "NS1" | CI.method == "NS2" ) {
    NC = (CI.method == "NS2")
    pXY = p.XY(x, y, n.x, n.y)
    if (LT==FALSE) {
      CI = polyroot2(HM.coef, AUC.hat=A, n.x=n.x, n.y=n.y, pXY=pXY, alpha=alpha, r=0, NC=NC)
    } else{
      if (start[2]==1) {
        ub = 1; start = start[1]
        lb = multiroot2(HMS.equation, start, AUC.hat=A, n.x=n.x, n.y=n.y, pXY=pXY, alpha=alpha, LT=LT, NC=NC, rtol = 1e-10, atol = 1e-10,...)$root  #original scale    
        CI = c(lb, ub)
      } else {
        CI = multiroot2(HMS.equation, start, AUC.hat=A, n.x=n.x, n.y=n.y, pXY=pXY, alpha=alpha, LT=LT, NC=NC, rtol = 1e-10, atol = 1e-10,...)$root  #original scale
      }      
    }
    CI.LT = logit(CI, LT=LT)
    if (variance) {V.LT = V.HM(A, n.x, n.y,  pXY=pXY, LT=LT, NC=NC,...)}
  }
  
  else if (CI.method == "NW") {
    V.LT = V.HM(A, n.x, n.y, NC=TRUE, LT=LT,...)
    if (is.na(V.LT)) {
      warning("Variance unestimable, the number of diseased or nondiseased subjects is less than or equal to one")
      CI.LT = c(NaN, NaN)
    } else {
      if (V.LT < 0) { 
        warning("Negative estimated variance set to zero")
        V.LT = 0
      }
      CI.LT = CI.base(A.LT,V.LT,alpha)
    }
  }
  
  else if (CI.method == "CM") {
    V.LT = V.CM(A, n.x, n.y, LT=LT,...)
    CI.LT = CI.base(A.LT,V.LT,alpha)
  }
  
  else if (CI.method == "RG") {
    probit = probit.RG(x,y,n.x,n.y)
    d = probit$delta; V = probit$V
    A = pnorm(d)
    d.CI = d + c(-1,1)*qnorm(1-alpha/2)*sqrt(V)
    CI = pnorm(d.CI)
    CI.LT = logit(CI, LT=LT) # logit transformation not available for RG (only for presentation purpose)
    V.LT = NA                # variance not available for RG(The V above is variance in probit scale)
  }
  
  else if (CI.method == "Bm") {
    V.LT = V.Bm(x, y, n.x, n.y, A, LT=LT,...)
    CI.LT = CI.base(A.LT,V.LT,alpha)
  }

  else if (CI.method == "Mee") {
    N.J.hat = Mee.stat(x, y, n.x, n.y)$N.J.hat
    if (LT==FALSE) {
      CI = polyroot2(Mee.coef, AUC.hat=A, N.J.hat=N.J.hat, alpha=alpha, r=0)
    } else{
      CI = multiroot2(Mee.equation, start, AUC.hat=A, N.J.hat=N.J.hat, alpha=alpha, LT=LT, rtol = 1e-10, atol = 1e-10,...)$root
    }
    CI.LT = logit(CI, LT=LT)
    if (variance) {
      V.LT=V.Mee(x, y, n.x=n.x, n.y=n.y, LT=LT,...) 
    }
  }
  
  else if (CI.method == "DB") {     #variance(or alp) has an unstable solution
    CI = DB(A, n.x, n.y, alpha, LT=LT,...)
    CI.LT = logit(CI, LT=LT)
    if (variance) {
      alp = multiroot2(DB.equation, start, AUC.hat=A, n.x=n.x, n.y=n.y, alpha=1,LT=LT)$root[1]
      V.LT = V.DB(alp, n.x, n.y, LT=LT)
    }
  }
  
  else if (CI.method == "DG") {  #variance(or alp) has an unstable solution
    CI = DG(A, n.x, n.y, alpha, LT=LT,...)
    CI.LT = logit(CI, LT=LT)
    if (variance) {
      delta = multiroot2(DG.equation, start, AUC.hat=A, n.x=n.x, n.y=n.y, alpha=1,LT=LT, rtol = 1e-10, atol = 1e-10,...)$root[1]
      V.LT = V.DG(delta, n.x, n.y, LT=LT,...)
    }
  }
  
  else if (CI.method == "DL") {
    V.LT = V.DeLong(x, y, n.x, n.y, A, LT=LT,...)
    CI.LT = CI.base(A.LT,V.LT,alpha)
  }
  
  else if (CI.method == "RJ"){  # logit (Rubin) with Jeffrey's prior
    Meestat <- Mee.stat(x,y)
    N.total <- Meestat$N.J.hat
    AUC.hat <- Meestat$AUC.hat
    N.success <- AUC.hat * N.total
    AUC.tilde <- (N.success + .5) / (N.total + 1)
    V.LT = 1 / ((N.total + 1) * AUC.tilde * (1 - AUC.tilde))
    CI.LT = CI.base(logit(AUC.tilde), V.LT, alpha)
    LT = TRUE      # fixing LT = TRUE only.
    A = AUC.tilde
  }

  ##### WilsonScore and ClopperPearson CI.methods are excluded from this study ##########
  else if (CI.method == "WS") {
    z.a2 = qnorm(1-alpha/2)
    CI = (2*n*A + z.a2^2 +cc*c(-1,1) +c(-1,1)* z.a2* sqrt(z.a2^2 + cc*c(-2,2) - 1/n*cc + 4*A*(n*(1-A) +1*cc)))/(2*(n+z.a2^2))
  }
  
  # V for CP???
  else if (CI.method == "CP") {
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
    CI.i$AUC.hat = fun(data, CI.method = CI.method[1],...)$AUC.hat
    for (i in 1:length(CI.method)) { CI.i[,i+1] = fun(data, CI.method = CI.method[i], ...)$CI }
  } else if (type == "landscape2") {
    CI.i = as.data.frame(matrix(NA,1,length(CI.method)*2+1))
    names(CI.i) = c("AUC.hat",paste0(rep(CI.method,each=2),c(".lb",".ub")))
    CI.i$AUC.hat = fun(data, CI.method = CI.method[1],...)$AUC.hat
    for (i in 1:length(CI.method)) { CI.i[1,(2*i):(2*i+1)] = fun(data, CI.method = CI.method[i], ...)$CI }  
  } else if (type=="portrait") {
    CI.i = as.data.frame(matrix(NA,length(CI.method)+1,2))
    names(CI.i) = c("lowerbound", "upperbound")
    rownames(CI.i)  = c("AUC.hat", CI.method)
    CI.i[1,] = fun(data, CI.method = CI.method[1],...)$AUC.hat
    for (i in 1:length(CI.method)) { CI.i[i+1,] = fun(data, CI.method = CI.method[i], ...)$CI }  
  } else {stop("Form is none of landscape1, landscape2, or portrait")}    
  return(CI.i)
}
