### Part I. CI functions
## 1-1b [CI] Auxillary functions: 1.1b.1 Rubin, 1.1.2 Error handling(mice, numerical solutions), 
library(mice)       # for Multiple Imputation
library(norm)         # for Multiple Imputation


# 1.1b.1 Combining Multiple imputations
# variance estimator for multiple imputation
Rubin = function(W, MI, alpha=0.05, print.r=FALSE, print.nu=FALSE) {
  if (is.na(MI[1])) {     # in case it is not Multiple imputation, bypass Rubin!
    v.final = W
    z.val=qnorm(1-alpha/2)
    r = 0
    nu = NA
  } else {
    m = length(MI)
    B = var(MI)
    v.final = W + B*(1+1/m)
    r = v.final / W  - 1   # missing ratio
    nu = (m-1)*(1+1/r)
    z.val = qt(1-alpha/2, nu)
  }
  result = data.frame(v.final=v.final, z.val=z.val)
  if (print.r) {result$r = r}
  if (print.nu) {result$nu = nu}
  return(result)
}


# Handling errors

# rootSolve:::multiroot without error
multiroot2 = function(f, start,..., silent=FALSE) {
  tmp <- try(multiroot(f, start,...),silent=silent)
  len <- length(start)
  if (class(tmp)=="try-error") {result <- list(); result$root <- rep(NA,len)} else {result <- tmp}
  return(result)
}

optim2 = function(equation, start, AUC.hat, method, lower, upper,..., silent=FALSE) {
  equation2 = function(...) {equation(...)^2}
  tmp <- try(mice(data = data, m = m, ...),silent=silent)
  lb <- try(optim(start[1], equation2, ..., AUC.hat=AUC.hat, lower=0, upper=AUC.hat, method="L-BFGS-B")$par, silent=silent)
  ub <- try(optim(start[2], equation2, ..., AUC.hat=AUC.hat, lower=AUC.hat, upper=1, method="L-BFGS-B")$par, silent=silent)
  if (class(lb)=="try-error") {lb <- NA}
  if (class(ub)=="try-error") {ub <- NA}
  return(c(lb,ub))
}

polyroot2 = function(fun,AUC.hat,..., silent=TRUE) {
  solutions = try(polyroot(fun(AUC.hat=AUC.hat,...)), silent=silent)
  if (class(solutions)=="try-error") {
    lb <- ub <- NA
  } else {
    sol.real = Re(solutions[round(Im(solutions),6)==0])
    sol.real = round(sol.real,6) # to deal with the problem where 1 is not actually 1 but 1-0.000000001
    lb = max(sol.real[sol.real< (AUC.hat)]) # unless AUC.hat==0, lb is strictly less than AUC.hat
    ub = min(sol.real[sol.real>=AUC.hat])   # ub is greater than AUC.hat with equality iff AUC.hat==1
    if (!is.finite(lb)) {lb=NA}; if (!is.finite(ub)) {ub=NA}
  }
  return(c(lb,ub))
}


# mice without error & complete
mice2 = function(data, m = 5, return.vec=c("diseaseR","marker"), ..., silent=FALSE) {
  data.imp <- try(mice(data = data, m = m, ...),silent=silent)
  if (class(data.imp) !="try-error") {
    data.comp = lapply(1:m, function(i) {complete(data.imp, action = i)[,return.vec]})
  } else {
    data.comp <- lapply(1:m, function(i) {NA})
  }
  return(data.comp)
}

# norm MI: adp.round is required for MI.norm
adp.round <- function(x, rounding="simple", vec=x) {
  # vec is the vector of a variable that is observed and imputed, and is only used for adaptive rounding.
  
  if (length(x) <= 1 & is.null(vec) & rounding == "adaptive") {warning("adaptive rounding is used only if the reference vector is provided")} 
  x = pmax(pmin(x,1),0)  # limiting the prediction range to [0,1]
  if (rounding == "simple") {
    x = round(x)
  } else if (rounding == "coinflip") {
    x = sapply(x, function(x) rbinom(1,1,x))
  } else if (rounding == "adaptive") {
    omega = mean(vec, na.rm = TRUE)
    cutoff = omega - qnorm(omega) * sqrt(omega*(1-omega))
    x = (x>=cutoff)*1
  } else {stop("rounding scheme should be either simple, coinflip, or adaptive.")}
  return(x)
}

MI.norm <- function(data, m=5, rounding="simple", rnd.vec = "diseaseR", return.vec=c("diseaseR","marker"),
                    rngseed=abs(round(rnorm(1)*10^5))+1, showits=FALSE) {
  # preliminary settings
  rngseed(rngseed)
  data.mat <- data.matrix(data)
  s <- prelim.norm(data.mat)
  mle <- em.norm(s, showits=showits)
  
  data.comp = lapply(1:m, function(a) {imp.norm(s = s, theta = mle, x = data.mat)})   # multiple imputations
  # data.comp = lapply(data.comp, function(data) {data[,rnd.vec] = pmax(pmin(data[,rnd.vec],1),0); return(data)}) # limiting the prediction range to [0,1]
  data.comp = lapply(data.comp, function(data) {data[,rnd.vec] = adp.round(data[,rnd.vec], rounding=rounding); return(data[,return.vec])}) # limiting the prediction range to [0,1]
}

