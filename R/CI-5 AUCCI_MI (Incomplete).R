### Part I. CI functions
## 1-5 [CI] AUCCI_MI(Incomplete):  1.5.1 AUCCI.MI

# Library ##############################################################################
# [package] mice, norm, [functoins] mice2, adp.norm, MI.norm are required: CI-1b.auxillary functions.R

# needs functions: AUC, data2xy, AUCCI, HMS.equation

# MI methods
MI.methods = data.frame(functions=c(rep("mice2",2),rep("MI.norm",3)), methods=c("pmm","logreg","simple","coinflip","adaptive"))
########################################################################################


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
      if (identical(MI.function, mice2)) {
        data.comp <- MI.function(data=data[,-1], method=MI.method, m=m, printFlag=FALSE,...)
      } else if (identical(MI.function, MI.norm)) {
        data.comp <- MI.function(data=data[,-1], rounding=MI.method, m=m, showits=FALSE,...)
      }
    }
    
    # 3. Otherwise, stop running the code.
    else {stop("data is neither a dataframe nor a list of imputed datasets.")}
  }
  
  ### B. Inference stage: point estimate
  {
    mi.stat = data.frame(theta = rep(NA, m), var = rep(NA, m))
    n.x.vec = rep(NA, m)
    for (i in 1:m) {
      temp = data2xy(data.comp[[i]],disease="diseaseR", marker="marker")
      x = temp$x ; y = temp$y
      mi.stat$theta[i] = AUC(x, y)
      n.x.vec[i] = length(x)
      if ("Mee" %in% CI.method) {
        mi.stat$N.J.hat[i] = Mee.stat(x, y)$N.J.hat
      }
    }
    n = dim(data.comp[[1]])[1]
    n.x = round(mean(n.x.vec, na.rm=TRUE))
      if (n.x <= 1 & mean(n.x.vec, na.rm=TRUE) >1) {n.x <- 2} # if n.x is 1, variance is not estimable.
    n.y = n - n.x
    mi.stat$theta.LT = logit(mi.stat$theta,LT=LT)
    mi.stat$div.LT = ifelse(LT,(mi.stat$theta[i]*(1-mi.stat$theta[i]))^2,1)
    mean.MI = mean(mi.stat$theta, na.rm=TRUE)
    mean.MI.LT = mean(mi.stat$theta.LT, na.rm=TRUE)
    var.B.LT = var(mi.stat$theta.LT, na.rm=TRUE)    # => Erase when not necessary any more after using funtion Rubin
    
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
      # when there is only one 0 or 1, V will be Inf (actually it's not inf, but nonestimable. So eliminate them when averaging.)
      mi.var.fin = mi.stat$var.LT[is.finite(mi.stat$var.LT)]
      # If every MI cases has unestimable variance(NA or Inf), Wald type CI can not be gotten.
      if (length(mi.var.fin) > 0) {
        Rubin = Rubin(W=mean(mi.var.fin), MI=mi.stat$theta.LT, alpha=alpha, print.nu = TRUE)
        var.MI.LT = Rubin$v.final         # W + (1/m)*B for MI
        if (is.na(var.MI.LT)) {
          CI <- CI.LT <- c(NA,NA)
        } else {
          nu.LT = Rubin$nu             # degree of freedom for MI
          if (is.na(nu.LT) & var.MI.LT==0) {nu.LT=Inf}
          CI.LT = CI.base(mean.MI.LT, var.MI.LT, alpha, qt, df=nu.LT)    # logit scale (will back-transform in the end)
          CI = expit(CI.LT, LT=LT)
        }
      } else {
        CI.LT <- CI <- c(NA,NA)
      }
    }
    
    else if (CI.method %in% CI.Score) {
      start = c( (mean.MI+.5)/2, (mean.MI+1)/2)
      start = mean.MI + c(-1,+1)*qnorm(1-alpha/2)*sqrt(mean.MI*(1-mean.MI)/n.x/n.y)
      if (score.MI == "fixed.r" & LT==FALSE) {  #if n.x==1, var is not estimable, make it population var to estimate r.
        W = V.HM(AUC=mean.MI, n.x=n.x, n.y=n.y, LT=LT, NC=FALSE, sample.var=ifelse(n.x==1,FALSE,TRUE))
        r = Rubin(W=W, MI=mi.stat$theta.LT, alpha=alpha, print.r=TRUE)$r
        mi.theta.fin = mi.stat$theta.LT[is.finite(mi.stat$theta.LT)]
        if (length(mi.theta.fin)==0) {
          CI.LT <- CI <- c(NA,NA)
        } else {
        if (all(mi.theta.fin==1)) {r=1}  #if all theta_(i)'s are 1, W=0, r=undefined, V=0
        if (CI.method == "NS1") {
          CI = polyroot2(HM.coef,AUC.hat=mean.MI, n.x=n.x, n.y=n.y, NC=FALSE, alpha=alpha, r=r)
        } else if (CI.method == "NS2") {
          CI = polyroot2(HM.coef,AUC.hat=mean.MI, n.x=n.x, n.y=n.y, NC=TRUE, alpha=alpha, r=r)
        } else if (CI.method == "Mee") {
          CI = polyroot2(Mee.coef,AUC.hat=mean.MI, N.J.hat=mean(mi.stat$N.J.hat, na.rm=TRUE), alpha=alpha, r=r)
        }
        }
      }
      else {
        if (CI.method == "NS1") {
          CI = multiroot2(HMS.equation, start, AUC.hat=mean.MI, n.x=n.x, n.y=n.y, alpha=alpha, MI=mi.stat$theta.LT, LT=LT, score.MI=score.MI, sample.var=FALSE, rtol = 1e-10, atol = 1e-10, ...)$root    
        } else if (CI.method == "NS2") {
          CI = multiroot2(HMS.equation,  start, AUC.hat=mean.MI, n.x=n.x, n.y=n.y, alpha=alpha, MI=mi.stat$theta.LT, LT=LT, NC=TRUE, score.MI=score.MI, sample.var=FALSE, rtol = 1e-10, atol = 1e-10, ...)$root 
        } else if (CI.method == "Mee") {
          CI = multiroot2(Mee.equation,  start, AUC.hat=mean.MI, N.J.hat = mean(mi.stat$N.J.hat, na.rm=TRUE), alpha=alpha, MI=mi.stat$theta.LT, LT=LT, score.MI=score.MI, sample.var=FALSE, rtol = 1e-10, atol = 1e-10, ...)$root
        }
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
            MI.DB$alp[i] = multiroot2(DB.equation, c(1,3), AUC.hat=mi.stat$theta[i], n.x=n.x, n.y=n.y, alpha=1 ,LT=LT, rtol = 1e-10, atol = 1e-10)$root[1]
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
            MI.DG$delta[i] = multiroot2(DG.equation, c(1,3), AUC.hat=mi.stat$theta[i], n.x=n.x, n.y=n.y, alpha=1,LT=LT, rtol = 1e-10, atol = 1e-10)$root[1]
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
