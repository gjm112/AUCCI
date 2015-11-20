### Part I. CI functions
## 1-1b [CI] Auxillary functions: 1.1b.1 Rubin, 1.1.2 Error handling(mice, numerical solutions), 

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
# mice without error
mice2 = function(data, m = 5, ..., silent=FALSE) {
  tmp <- try(mice(data = data, m = m, ...),silent=silent)
  if (class(tmp)=="try-error") {
    result <- "error"} else {result <- tmp}
  return(result)
}

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