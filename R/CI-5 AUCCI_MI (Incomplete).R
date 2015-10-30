### Part I. CI functions
## 1-5 [CI] AUCCI_MI(Incomplete):  1.5.1 AUCCI.MI

# Library ##############################################################################
library(mice)       # for Multiple Imputation
library(norm)         # for Multiple Imputation

# needs functions: AUC, data2xy, AUCCI, HMS.equation

# MI methods
MI.methods = data.frame(functions=c(rep("mice",2),rep("norm",1)), methods=c("pmm","logreg","imp.norm"))
########################################################################################

MI.norm <- function(data) {
  data.mat <- data.matrix(data)
  s <- prelim.norm(data.mat)
  mle <- em.norm(s)
  imp.norm(s = s, theta = mle, x = data.mat)
}

## 1.5.1 AUCCI.MI ##################################################################### 
# For score types DO NOT use LT!!! very unstable for theta near 1

AUCCI.MI = function(data, MI.function, MI.method, score.MI = "fixed.r", m, CI.method, alpha=.05, LT=FALSE, variance=FALSE, MI.thetas=FALSE, ...) {
  # note: predictorMatrix(for mice) should exclude disease
  # LT: logit transformation      # not available for RG
  # returning variance when logit=T is var(logit(theta))
  # score.MI = type of adding B to W: "fixed.r", "fixed.B"
  
  ### A. imputation stage : output = data.comp and m
  {
    # 1. when multiple imputation is already done and stored in a dataset as a list.
    # For efficient simulation
    if (class(data) == "list") {
      m = length(data)
      data.comp = data
    }
    
    # 2. when multiple imutation needs to be done for a dataframe
    # For convenience of one-time CI construction. Very unefficient for large size of simulation
    else if (class(data) == "data.frame") {
      if (identical(MI.function, mice)) {
        data.imp <- MI.function(data=data, method=MI.method, m=m, predictorMatrix=cbind(0,(1 - diag(1, ncol(data)))[,-1]), printFlag=FALSE,...)
        data.comp <- list()
        for (i in 1:m) {data.comp[[i]] = complete(data.imp, action = i)[,c("diseaseR", "marker")]}
      } else if (identical(MI.function, norm)) {
        #######  "TBD"  
        
      }
    }
    
    # 3. Otherwise, stop running the code.
    else {stop("data is neither a dataframe nor a list of imputed datasets.")}
  }
  
  ### B. Inference stage: point estimate
  {
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
    
  }
  
  ### C. CI stage
  {
    # Categorizing CI types
    CI.Wald = c("HM1", "HM2", "NW", "CM", "Bm", "DL")
    CI.Score = c("NS1", "NS2", "Mee")
    CI.root = c("DB", "DG")
    CI.other = c("RG")
    
    if (CI.method %in% CI.Wald) {
      for (i in 1:m) {
        mi.stat$var.LT[i] = AUCCI(data.comp[[i]], CI.method=CI.method, disease="diseaseR", alpha=alpha, variance = TRUE, LT=LT, ...)$V.hat
      }
      Rubin = Rubin(W=mean(mi.stat$var.LT), MI=mi.stat$theta.LT, alpha=alpha, print.nu = TRUE)
      var.MI.LT = Rubin$v.final         # W + (1/m)*B for MI
      nu.LT = Rubin$nu             # degree of freedom for MI
      CI.LT = CI.base(mean.MI.LT, var.MI.LT, alpha, qt, df=nu.LT)    # logit scale (will back-transform in the end)
      CI = expit(CI.LT, LT=LT)
    }
    
    else if (CI.method %in% CI.Score) {
      start = c( max(mean.MI - 0.1, mean.MI/2), min(mean.MI + 0.1, (mean.MI+1)/2))
      if (CI.method == "NS1") {
        CI = multiroot(HMS.equation, start, AUC.hat=mean.MI, n.x=n.x, n.y=n.y, alpha=alpha, MI=mi.stat$theta.LT, LT=LT, score.MI=score.MI, sample.var=FALSE, rtol = 1e-10, atol = 1e-10, ...)$root}
      else if (CI.method == "NS2"){
        CI = multiroot(NC.equation,  start, AUC.hat=mean.MI, n.x=n.x, n.y=n.y, alpha=alpha, MI=mi.stat$theta.LT, LT=LT, score.MI=score.MI, sample.var=FALSE, rtol = 1e-10, atol = 1e-10, ...)$root}  
      else if (CI.method == "Mee") {
        CI = multiroot(Mee.equation,  start, AUC.hat=mean.MI, x=x, y=y, n.x=n.x, n.y=n.y, alpha=alpha, MI=mi.stat$theta.LT, LT=LT, score.MI=score.MI, sample.var=FALSE, rtol = 1e-10, atol = 1e-10, ...)$root
      }
      CI.LT = logit(CI,LT=LT)
      if (variance) {
        for (i in 1:m) {
          mi.stat$var.LT[i] = AUCCI(data.comp[[i]], CI.method=CI.method, disease="diseaseR", alpha=alpha, variance = TRUE, LT=LT, ...)$V.hat
        }
        Rubin = Rubin(W=mean(mi.stat$var.LT), MI=mi.stat$theta.LT, alpha=alpha)
        var.MI.LT = Rubin$v.final         # W + (1/m)*B for MI
      }
    }
    
    else if (CI.method %in% CI.root) {
      if (CI.method == "DB") {
        CI = DB(AUC.hat=mean.MI, n.x=n.x, n.y=n.y, alpha=alpha, MI=mi.stat$theta.LT, LT=LT)
        CI.LT = logit(CI, LT=LT)
        if (variance) {
          MI.DB = data.frame(alp = rep(NA, m), var = rep(NA, m))
          for (i in 1:m) {
            MI.DB$alp[i] = multiroot(DB.equation, c(1,3), AUC.hat=mi.stat$theta[i], n.x=n.x, n.y=n.y, alpha=1 ,LT=LT, rtol = 1e-10, atol = 1e-10)$root[1]
            MI.DB$var.LT[i] = V.DB(MI.DB$alp[i], n.x, n.y, LT=LT, ...)
          }
          Rubin = Rubin(W=mean(MI.DB$var.LT), MI=mi.stat$theta.LT, alpha=alpha)
          var.MI.LT = Rubin$v.final         # W + (1/m)*B for MI
        }
      }
      
      else if (CI.method == "DG") {
        CI = DG(AUC.hat=mean.MI, n.x=n.x, n.y=n.y, alpha=alpha, MI=mi.stat$theta.LT, LT=LT, ...)
        CI.LT = logit(CI, LT=LT)
        if (variance) {
          MI.DG = data.frame(delta = rep(NA, m), var.LT = rep(NA, m))
          for (i in 1:m) {
            MI.DG$delta[i] = multiroot(DG.equation, c(1,3), AUC.hat=mi.stat$theta[i], n.x=n.x, n.y=n.y, alpha=1,LT=LT, rtol = 1e-10, atol = 1e-10)$root[1]
            MI.DG$var.LT[i] = V.DG(MI.DG$delta[i], n.x, n.y, LT=LT, ...)
          }
          Rubin = Rubin(W=mean(mi.stat$var.LT), MI=mi.stat$theta.LT, alpha=alpha)
          var.MI.LT = Rubin$v.final         # W + (1/m)*B for MI
        }
      }
    }  
    else if (CI.method == "RG") {         # LT, variance not available
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
    
    else {warning(paste(CI.method,"is not available"))} 
  }
  
  result = list(AUC.hat = mean.MI, CI = CI);  if (variance) {result$V.hat = var.MI.LT}; if (MI.thetas) {result$theta.MI = mi.stat$theta; result$theta.LT = mi.stat$theta.LT }
  return(result)
}