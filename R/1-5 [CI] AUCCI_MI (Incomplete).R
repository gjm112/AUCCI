### Part I. CI functions
## 1-5 [CI] AUCCI_MI(Incomplete):  1.5.1 AUCCI.MI

# Library ##############################################################################
library(mice)       # for Multiple Imputation
library(norm)         # for Multiple Imputation

# needs functions: AUC, data2xy, AUCCI, HMS.equation

# MI methods
MI.methods = data.frame(functions=c(rep("mice",2),rep("norm",1)), methods=c("pmm","logreg","imp.norm"))
########################################################################################


## 1.5.1 AUCCI.MI ##################################################################### 

## preset arguments
data = temp
m = 10  # number of imputations
MI.function = mice
MI.method = "pmm"
alpha = .05

AUCCI.MI = function(data, MI.function, MI.method, m, CI.method, alpha=.05, LT=FALSE, variance=FALSE, MI.thetas=FALSE, ...) {
  # note: predictorMatrix(for mice) should exclude disease
  # LT: logit transformation      # not available for RG
  # A. imputation stage
  # returning variance when logit=T is var(logit(theta))
  if (identical(MI.function, mice)) {
    data.imp <- MI.function(data=data, method=MI.method, m=m, predictorMatrix=cbind(0,(1 - diag(1, ncol(data)))[,-1]), printFlag=FALSE,...)
    data.comp <- list()
    for (i in 1:m) {data.comp[[i]] = complete(data.imp, action = i)[,c("diseaseR", "marker")]}
  } else if (identical(MI.function, norm)) {
    #######  "TBD"  
    
  }
  
  # B. Inference stage: point estimate
  mi.stat = data.frame(theta = rep(NA, m), var = rep(NA, m))
  for (i in 1:m) {
    temp = data2xy(data.comp[[i]],disease="diseaseR", marker="marker")
    x = temp$x ; y = temp$y
    mi.stat$theta[i] = AUC(x, y)
  }
  n.x = length(x); n.y = length(y); n = n.x + n.y
  mi.stat$theta.LT = logit(mi.stat$theta,LT=LT)
  mi.stat$div.LT = ifelse(LT,(mi.stat$theta[i]*(1-mi.stat$theta[i]))^2,1)
  mean.MI = mean(mi.stat$theta)
  mean.MI.LT = mean(mi.stat$theta.LT)
  var.B.LT = var(mi.stat$theta.LT)    # => Erase when not necessary any more after using funtion Rubin
  
  # C. CI stage
  # Categorizing CI types
  CI.Wald = c("HanleyMcNeilWald", "HanleyMcNeilExponential", "NewcombeExponential", 
              "CortesMohri", "Bamber", "MannWhitney", "DeLong")
  CI.Score = c("HanleyMcNeilScore", "NewcombeScore", "HalperinMee", "WilsonScore")
  CI.root = c("DoubleBeta", "DoubleGaussian")
  CI.Gauss = c("ReiserGuttman")
  CI.other = c("MannWhitneyLT", "ClopperPearson")
  
  if (CI.method %in% CI.Wald) {
    for (i in 1:m) {
      mi.stat$var.LT[i] = AUCCI(data.comp[[i]], method=CI.method, disease="diseaseR", alpha=alpha, variance = TRUE, LT=LT)$V.hat
    }
    Rubin = Rubin(W=mean(mi.stat$var.LT), MI=mi.stat$theta.LT, alpha=alpha, print.nu = TRUE)
    var.MI.LT = Rubin$v.final         # W + (1/m)*B for MI
    nu.LT = Rubin$nu             # degree of freedom for MI
    CI.LT = CI.base(mean.MI.LT, var.MI.LT, alpha, qt, df=nu.LT)    # logit scale (will back-transform in the end)
  }
  if (CI.method %in% CI.Score) {
    if (CI.method == "HanleyMcNeilScore") {
      CI = multiroot(HMS.equation, c(0.5,0.9), AUC.hat=mean.MI, n.x=n.x, n.y=n.y, alpha=alpha, MI=mi.stat$theta.LT, LT=LT)$root}
    else if (CI.method == "NewcombeScore") {
      CI = multiroot(NC.equation, c(0.5,0.9), AUC.hat=mean.MI, n.x=n.x, n.y=n.y, alpha=alpha, MI=mi.stat$theta.LT, LT=LT)$root}  
    else if (CI.method == "HalperinMee") {
      CI = multiroot(Halperin.equation, c(0.5,0.9), AUC.hat=mean.MI, x=x, y=y, n.x=n.x, n.y=n.y, alpha=alpha, MI=mi.stat$theta.LT, LT=LT)$root
    }
    
    CI.LT = logit(CI,LT=LT)
    if (variance) {
      for (i in 1:m) {
        mi.stat$var.LT[i] = AUCCI(data.comp[[i]], method=CI.method, disease="diseaseR", alpha=alpha, variance = TRUE, LT=LT)$V.hat
      }
      Rubin = Rubin(W=mean(mi.stat$var.LT), MI=mi.stat$theta.LT, alpha=alpha)
      var.MI.LT = Rubin$v.final         # W + (1/m)*B for MI
    }
  }
  if (CI.method %in% CI.root) {
    if (CI.method == "DoubleBeta") {
      CI = DB(AUC.hat=mean.MI, n.x=n.x, n.y=n.y, alpha=alpha, MI=mi.stat$theta.LT, LT=LT)
      CI.LT = logit(CI, LT=LT)
      if (variance) {
        MI.DB = data.frame(alp = rep(NA, m), var = rep(NA, m))
        for (i in 1:m) {
          MI.DB$alp[i] = multiroot(DB.equation, c(1,3), AUC.hat=mi.stat$theta[i], n.x=n.x, n.y=n.y, alpha=1 ,LT=LT)$root[1]
          MI.DB$var.LT[i] = V.DB(MI.DB$alp[i], n.x, n.y, LT=LT)
        }
        Rubin = Rubin(W=mean(MI.DB$var.LT), MI=mi.stat$theta.LT, alpha=alpha)
        var.MI.LT = Rubin$v.final         # W + (1/m)*B for MI
      }
    }
    else if (CI.method == "DoubleGaussian") {
      CI = DG(AUC.hat=mean.MI, n.x=n.x, n.y=n.y, alpha=alpha, MI=mi.stat$theta.LT, LT=LT)
      CI.LT = logit(CI, LT=LT)
      if (variance) {
        MI.DG = data.frame(delta = rep(NA, m), var.LT = rep(NA, m))
        for (i in 1:m) {
          MI.DG$delta[i] = multiroot(DG.equation, c(1,3), AUC.hat=mi.stat$theta[i], n.x=n.x, n.y=n.y, alpha=1,LT=LT)$root[1]
          MI.DG$var.LT[i] = V.DG(MI.DG$delta[i], n.x, n.y, LT=LT)
        }
        Rubin = Rubin(W=mean(mi.stat$var.LT), MI=mi.stat$theta.LT, alpha=alpha)
        var.MI.LT = Rubin$v.final         # W + (1/m)*B for MI
      }
    }
  }  
  if (CI.method == "ReiserGuttman") {         # LT, variance not available
    MI.RG = data.frame(delta = rep(NA, m), var = rep(NA, m))
    for (i in 1:m) {
      temp = data2xy(data.comp[[i]],disease="diseaseR", marker="marker")
      x = temp$x ; y = temp$y
      MI.RG$delta[i] = probit.RG(x, y)$delta
      MI.RG$var[i] = probit.RG(x, y)$V
    }
    delta.mean = mean(MI.RG$delta)
    Rubin = Rubin(W=mean(MI.RG$var), MI=MI.RG$delta, alpha=alpha)
    var.MI = Rubin$v.final         # W + (1/m)*B for MI
    z.val = Rubin$z.val
    d.CI = delta.mean + c(-1,1)*z.val*sqrt(var.MI)
    CI = pnorm(d.CI)
    CI.LT = logit(CI, LT=LT) # logit transformation not available for RG (only for presentation purpose)
    V.MI.LT = NA                # variance not available for RG(The V above is variance in probit scale)
    mean.MI=pnorm(delta.mean)
    var.MI.LT = NA
  }
  
  result = list(AUC.hat = mean.MI, CI = CI);  if (variance) {result$V.hat = var.MI.LT}; if (MI.thetas) {result$theta.MI = mi.stat$theta; result$theta.LT = mi.stat$theta.LT }
  return(result)
}

# Example
AUCCI.MI(data[1:100,], MI.function=mice, MI.method="pmm", alpha=.05, CI.method="HanleyMcNeilWald", m = m)
CI.i(data[1:100,], fun=AUCCI.MI, CI.method=CI.methods, type="landscape2",MI.function=mice, MI.method="pmm", alpha=.05, m = m)

# AUCCI(data,"HanleyMcNeilWald")
