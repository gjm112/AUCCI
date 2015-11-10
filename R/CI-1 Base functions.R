### Part I. CI functions
## 1-1 [CI] Base functions: 1.1.1 CI.base, 1.1.2 Variance estimators

# Library ##############################################################################
library(rootSolve)  # for Newton-Raphson Method (0.2.4)
library(sn)         # for Owen's T-function in Double Gaussian Method (0.2.15)

# CI methods :    ref to 2-3.[SIM] Simulator, 3.evaluation
CI.methods = c("Bm", "HM1", "HM2", "NS1", "NW", "NS2", "Mee", "DL", "RG", "CM", "DB", "DG")
### Descriptions
# Bm : Bamber
# HM1 : Hanley-McNeil 1 (Wald, nonparametric var)
# HM2 : Hanley-McNeil 2 (Wald, parametric var:exponential)
# NS1 : Newcombe Score1 (Score + HM2)
# NW : Newcombe Wald (HM2 + balanced design)
# NS2 : Newcombe Score 2 (NS1 + balanced design)
# Mee : Mee (Score)
# DL : DeLong (Wald)
# RG : Reiser-Guttman (Parametric: Gaussian)
# DB : Double beta (Parametric: double beta)
# DG : Double Gaussian (Parametric: double Gaussian)
# CM : Cortes-Mohri (Wald, Combinatoric variance)

###               "CP","WS" => excluded from this study: not AUC specific
########################################################################################

## 1.1.1 Basic internal functions ######################################################
# CI.base, AUC.classic, AUC
CI.base <- function(mu.hat, Var, alpha, dist=qnorm,...) {
  return(mu.hat + c(-1,1) * dist(1-alpha/2,...)*sqrt(Var))
}
AUC.classic <- function(x, y, n.x=length(x), n.y=length(y)) {
  A <- 0; for (i in 1:n.x) { A = A + ( sum(y > x[i]) + sum(y==x[i])/2 ) }
  return(A = A/ n.x / n.y)
}
# new AUC function added with much improved efficiency
AUC <- function(x, y, n.x=length(x), n.y=length(y), data=NA, disease="disease", marker="marker") {
  if (!is.na(data[1])[1]) {
    xy = data2xy(data, disease=disease, marker=marker)
    x = xy$x
    y = xy$y
  }
  predictions = c(x,y)
  labels = c(rep(0,n.x), rep(1,n.y))
  pred.order = order(predictions, decreasing = TRUE)
  predictions.sorted = predictions[pred.order]
  tp = cumsum(labels[pred.order] == 1)
  fp = cumsum(labels[pred.order] == 0)
  dups <- rev(duplicated(rev(predictions.sorted)))
  tp <- c(0, tp[!dups])
  fp <- c(0, fp[!dups])
  #cutoffs <- c(Inf, predictions.sorted[!dups])
  # excerpted until here from <ROCR> prediction function
  tpr <- tp[-1]/max(tp)
  d.fp <- fp[-1] - fp[-length(fp)]
  d.fpr <- d.fp/max(fp)
  return(sum(tpr*d.fpr))
}

# data converter: from dataframe to x, y vectors
data2xy <- function(data, disease="disease", marker="marker") {
  x = data[data[,disease]==0, marker]
  y = data[data[,disease]==1, marker]
  return(list(x=x,y=y))
}

# logit, expit
logit <- function(data, LT=TRUE) {
  if (LT) {return(log(data/(1-data)))} else {return(data)}
}
expit <- function(data, LT=TRUE) {
  if (LT) {return((1+exp(data)^-1)^-1)} else {return(data)}
}

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

# rootSolve:::multiroot without error
multiroot2 = function(f, start,..., silent=FALSE) {
  tmp <- try(multiroot(f, start,...),silent=silent)
  len <- length(start)
  if (class(tmp)=="try-error") {result <- list(); result$root <- rep(NA,len)} else {result <- tmp}
  return(result)
}


## 1.1.2 variance estimators ###########################################################
# Q.stat for HM and its derivatives
Q.stat <- function(x, y, n.x=length(x), n.y=length(y)) {
  Q2 <- Q1 <- 0
  for (i in 1:n.x) {Q1 = Q1 + ( sum(y > x[i]))^2}; Q1 = Q1 /(n.x * n.y^2)  # no ties
  for (j in 1:n.y) {Q2 = Q2 + ( sum(y[j] > x))^2}; Q2 = Q2 /(n.x^2 * n.y)
  return(data.frame(Q1 = Q1, Q2 = Q2))
}
V.HM <- function(AUC, n.x, n.y, Q1 = AUC/(2-AUC), Q2=2* AUC^2 /(AUC + 1), LT=FALSE, NC=FALSE, sample.var=TRUE, ...) {
  # default: HME method
  N = (n.x+n.y)/2
  if (LT) {div.LT = (AUC*(1-AUC))^2} else {div.LT=1}
  nxny = ifelse(sample.var, (n.x-1)*(n.y-1), n.x*n.y)  # sample.var: dividing by (nx-1)(ny-1) is unbiased
  if (NC==FALSE) {V = (AUC * (1-AUC) + (n.y -1)*(Q1 - AUC^2) + (n.x-1)*(Q2 - AUC^2)) /nxny /div.LT}
  else {V = AUC*(1-AUC)/nxny * ( (2*N-1) - (3*N-3) / ((2-AUC)*(AUC+1)) )/div.LT}
  return(V)
}

### for score method it should be populatoin variance not sample variance (sample.var=FALSE).
HMS.equation = function(AUC, AUC.hat, n.x, n.y, alpha, MI=NA, LT=FALSE, NC=FALSE, sample.var=FALSE, score.MI = "fixed.r") {
  W = V.HM(AUC, n.x, n.y, LT=LT, NC=NC, sample.var=sample.var)
  if (score.MI == "fixed.r") {
    ## 1. fixed r
    W.hat = V.HM(AUC.hat, n.x, n.y, LT=LT, NC=NC, sample.var=sample.var)
    Rubin = Rubin(W=W.hat, MI=MI, alpha=alpha, print.r=TRUE)
    z.val = Rubin$z.val       # nu is fixed (since r is fixed)
    r = Rubin$r
    return((logit(AUC,LT=LT) - min(10,logit(AUC.hat,LT=LT)))^2 - z.val^2*W*(1+r) )  
  } else if (score.MI == "fixed.B") {
    ## 2. fixed B hat
    Rubin = Rubin(W=W, MI=MI, alpha=alpha)
    V = Rubin$v.final         # W + (1/m)*B for MI
    z.val = Rubin$z.val       # t value for MI
    return((logit(AUC,LT=LT) - min(10,logit(AUC.hat,LT=LT)))^2 - z.val^2*V)
  }
}
V.CM <- function(AUC, n.x, n.y, LT=FALSE){
  z1 = function (w, n.x, n.y, k) {1-0.5*(w/n.y+(k-w)/n.x)}
  z2 = function (w, n.x, n.y, k) 
    (n.x*w^2+n.y*(k-w)^2+n.x*(n.x+1)*w+n.y*(n.y+1)*(k-w)-2*w*(k-w)*(n.x+n.y+1))/(12*n.x^2*n.y^2)
  F.base <- function (w, n.x, n.y, k) {choose(n.x-k+2*w,w)*choose(n.y+k-2*w,k-w)}
  k=round((n.x+n.y)*(1-AUC))
  F1=F2=F3=F.den=0
  for (w in 0:k){
    F1 = F1 + F.base(w, n.x, n.y, k)*z1(w, n.x, n.y, k)^2
    F.den = F.den + F.base(w, n.x, n.y, k)
    F2 = F2 + F.base(w, n.x, n.y, k)*z1(w, n.x, n.y, k)
    F3 = F3 + F.base(w, n.x, n.y, k)*z2(w, n.x, n.y, k)
  }
  if (LT) {div.LT = (AUC*(1-AUC))^2} else {div.LT=1}
  return(V = (F1/F.den - F2^2/F.den^2 + F3/F.den)/div.LT)
}
probit.RG <- function(x, y, n.x=length(x), n.y=length(y)){
  mu.x = mean(x); mu.y = mean(y)
  sig.x = sd(x);  sig.y = sd(y); s2 = sig.x^2 + sig.y^2
  d = (mu.y-mu.x)/sqrt(s2)
  f = s2^2/((sig.x^4/(n.x-1))+(sig.y^4/(n.y-1)))
  M = s2/((sig.x^2/n.x)+(sig.y^2/n.y))
  V = 1/M + d^2/(2*f)
  return(data.frame(delta=d, V=V))
}

V.Bm <- function(x, y, n.x=length(x), n.y=length(y), AUC.hat=AUC(x,y,n.x,n.y), LT=FALSE, sample.var=TRUE){
  b.yyx <- b.xxy <- p.xy <- 0
  for (i in 1:n.x) {b.yyx = b.yyx + sum(y < x[i])^2 + sum(x[i] < y)^2 -2*sum(y < x[i])*sum(x[i] < y) }
  b.yyx = b.yyx/(n.x*n.y^2)
  for (j in 1:n.y) {b.xxy = b.xxy + sum(x < y[j])^2 + sum(y[j] < x)^2 -2*sum(x < y[j])*sum(y[j] < x) }
  b.xxy = b.xxy/(n.y*n.x^2)
  for (i in 1:n.x) {  p.xy = p.xy + sum(y != x[i] ) }
  p.xy = p.xy /(n.x*n.y)
  nxny = ifelse(sample.var, (n.x-1)*(n.y-1), n.x*n.y)  # sample.var: dividing by (nx-1)(ny-1) is unbiased
  if (LT) {div.LT = (AUC.hat*(1-AUC.hat))^2} else {div.LT=1}
  return(1/(4*nxny)*(p.xy + (n.y -1)*b.xxy + (n.x -1)*b.yyx - 4*(n.y + n.x -1)*(AUC.hat -0.5)^2)/div.LT)
}
# p.stat for Halperin Mee method
p.stat <- function(x, y, n.x=length(x), n.y=length(y)) {
  Q = Q.stat(x, y, n.x, n.y); Q1 <- Q$Q1;  Q2 <- Q$Q2
  y.mat <- matrix(rep(y,n.x),ncol=n.x)
  x.mat <- t(matrix(rep(x,n.y),nrow=n.x))
  U.sq <- sum(((y.mat > x.mat)+(y.mat == x.mat)/2)^2)
  p2 <- (n.x*n.y^2*Q1 - U.sq)/(n.x*n.y*(n.y-1))
  p1 <- (n.y*n.x^2*Q2 - U.sq)/(n.y*n.x*(n.x-1))
  return(data.frame(p1=p1,p2=p2))
}
##### after validation of new one, replace old one.
Mee.stat.old <- function(x, y, n.x=length(x), n.y=length(y)) {
  ranks = rank(c(x,y))
  rank.x <- ranks[1:n.x]
  rank.y <- ranks[-1:-n.x]
  AUC.hat <- AUC(x, y, n.x, n.y)
  increment = ifelse(AUC.hat < 0.5, +0.5, -0.5)
  while (min(AUC.hat,1-AUC.hat)*sqrt(n.x*n.y) < 0.5){
    rank.y = rank.y + increment
    x <- rank.x ; y <- rank.y
    AUC.hat = AUC(x, y, n.x, n.y)
  }
  p = p.stat(x, y, n.x, n.y)
  rho.hat = (p-AUC.hat^2)/(AUC.hat-AUC.hat^2)
  N.J.hat = n.x*n.y/(((n.x-1)*rho.hat$p1 +1)/(1-1/n.y) + ((n.y-1)*rho.hat$p2 + 1 )/(1-1/n.x) )
  return(data.frame(p1=p$p1, p2=p$p2, AUC.hat=AUC.hat, N.J.hat=N.J.hat))
}
### new one. replace the old one after validation!!
Mee.stat <- function(x, y, n.x=length(x), n.y=length(y)) {
  AUC.hat.0 <- AUC.hat <- AUC(x, y, n.x, n.y)
  increment = ifelse(AUC.hat < 0.5, +0.5, -0.5)
  while (min(AUC.hat,1-AUC.hat)*sqrt(n.x*n.y) < 0.5){
    y = y + increment
    AUC.hat = AUC(x, y, n.x, n.y)
  }
  p = p.stat(x, y, n.x, n.y)
  rho.hat = (p-AUC.hat^2)/(AUC.hat-AUC.hat^2)
  N.J.hat = n.x*n.y/(((n.x-1)*rho.hat$p1 +1)/(1-1/n.y) + ((n.y-1)*rho.hat$p2 + 1 )/(1-1/n.x) )
  return(data.frame(AUC.hat=AUC.hat.0, p1=p$p1, p2=p$p2, N.J.hat=N.J.hat))
}
V.Mee.old <- function(x, y, n.x=length(x), n.y=length(y), AUC=AUC(x,y,n.x,n.y), LT=FALSE, N.J.hat=NA, sample.var=TRUE) {
  # sample.var is not relevant. for compatibility purpose. (V.Mee is upward biased)
  if (is.na(N.J.hat)) {N.J.hat = Mee.stat(x,y,n.x,n.y)$N.J.hat}
  V = AUC*(1-AUC) / N.J.hat
  if (LT) {div.LT = (AUC*(1-AUC))^2} else {div.LT=1}
  return( V/div.LT )
}
Mee.equation.old <- function(AUC, AUC.hat, x, y, n.x, n.y, alpha, MI=NA, LT=FALSE, sample.var=TRUE, score.MI = "fixed.r"){
  W = V.Mee(x, y, n.x, n.y, AUC, LT=LT, sample.var=sample.var)
  if (score.MI == "fixed.r") {
    ## 1. fixed r
    W.hat = V.Mee(x, y, n.x, n.y, AUC.hat, LT=LT, sample.var=sample.var)
    Rubin = Rubin(W=W.hat, MI=MI, alpha=alpha, print.r=TRUE)
    z.val = Rubin$z.val       # nu is fixed (since r is fixed)
    r = Rubin$r
    return((logit(AUC,LT=LT) - logit(AUC.hat,LT=LT))^2 - z.val^2*W*(1+r) )  
  } else if (score.MI == "fixed.B") {
    Rubin = Rubin(W=W, MI=MI, alpha=alpha)
    V = Rubin$v.final         # W + (1/m)*B for MI
    z.val = Rubin$z.val       # t value for MI
    return((logit(AUC,LT=LT) - logit(AUC.hat,LT=LT))^2 - z.val^2*V)
  }
}

V.Mee <- function(AUC, N.J.hat, LT=FALSE, sample.var=TRUE) {
  # sample.var is not relevant. for compatibility purpose. (V.Mee is upward biased)
  V = AUC*(1-AUC) / N.J.hat
  if (LT) {div.LT = (AUC*(1-AUC))^2} else {div.LT=1}
  return( V/div.LT )
}
Mee.equation <- function(AUC, AUC.hat, N.J.hat, alpha, MI=NA, LT=FALSE, sample.var=TRUE, score.MI = "fixed.r"){
  W = V.Mee(AUC, N.J.hat, LT=LT, sample.var=sample.var)
  if (score.MI == "fixed.r") {
    ## 1. fixed r
    W.hat = V.Mee(AUC.hat, N.J.hat, LT=LT, sample.var=sample.var)
    Rubin = Rubin(W=W.hat, MI=MI, alpha=alpha, print.r=TRUE)
    z.val = Rubin$z.val       # nu is fixed (since r is fixed)
    r = Rubin$r
    return((logit(AUC,LT=LT) - logit(AUC.hat,LT=LT))^2 - z.val^2*W*(1+r) )  
  } else if (score.MI == "fixed.B") {
    Rubin = Rubin(W=W, N.J.hat, MI=MI, alpha=alpha)
    V = Rubin$v.final         # W + (1/m)*B for MI
    z.val = Rubin$z.val       # t value for MI
    return((logit(AUC,LT=LT) - logit(AUC.hat,LT=LT))^2 - z.val^2*V)
  }
}

V.DB <- function(alp, n.x, n.y, LT=FALSE) {
  R1 = gamma(alp + 1)^2 / gamma(2*alp + 1)
  R2 = gamma(2*alp +1)*gamma(alp +1)/gamma(3*alp+1)
  theta = 1 - R1
  Q = 1 - 2*R1 + R2
  if (LT) {div.LT = (theta*(1-theta))^2} else {div.LT=1}
  return(V = (theta*(1-theta) + (n.x+n.y-2)*(Q-theta^2))/n.x/n.y/div.LT)
}
DB.equation <- function(alp, AUC.hat, n.x, n.y, alpha, MI=NA, LT=FALSE){
  AUC = 1 - gamma(alp + 1)^2 / gamma(2*alp + 1)
  W = V.DB(alp, n.x, n.y, LT=LT)
  Rubin = Rubin(W=W, MI=MI, alpha=alpha)
  V = Rubin$v.final         # W + (1/m)*B for MI
  z.val = Rubin$z.val       # t value for MI
  return(logit(AUC,LT=LT) + c(+1,-1)*z.val*sqrt(V) - logit(AUC.hat,LT=LT))
}
DB <- function(AUC.hat, n.x, n.y, alpha,...) {
  start.1 <- min(3, max(-log(1-AUC.hat,2) - 1, 1))
  start.2 <- min(7, -log(1-AUC.hat,2) + 1)
  start <- c(start.1, start.2)
  alp <- multiroot2(DB.equation, start, AUC.hat=AUC.hat, n.x=n.x, n.y=n.y, alpha=alpha,...)$root
  return(1 - gamma(alp + 1)^2 / gamma(2*alp + 1))
}

V.DG <- function(delta, n.x, n.y, LT=FALSE){
  theta = pnorm(delta/sqrt(2))
  Owen = T.Owen(delta/sqrt(2), 1/sqrt(3))
  if (LT) {div.LT = (theta*(1-theta))^2} else {div.LT=1}
  return(V = ((n.x + n.y - 1)*theta*(1-theta) - 2*(n.x + n.y -2) * Owen ) / (n.x * n.y) /div.LT)
}
DG.equation <- function(delta, AUC.hat, n.x, n.y, alpha, MI=NA, LT=FALSE) {
  W = V.DG(delta, n.x, n.y, LT=LT)
  Rubin = Rubin(W=W, MI=MI, alpha=alpha)
  V = Rubin$v.final         # W + (1/m)*B for MI
  z.val = Rubin$z.val       # t value for MI
  AUC = pnorm(delta/sqrt(2))
  return(logit(AUC,LT=LT) + c(+1,-1)*z.val*sqrt(V) - logit(AUC.hat,LT=LT))
}
DG <- function(AUC.hat, n.x, n.y, alpha,...) {
  delta <- multiroot2(DG.equation, c(0.5, 0.9), AUC.hat=AUC.hat, n.x=n.x, n.y=n.y, alpha=alpha,...)$root
  return(pnorm(delta/sqrt(2)))
}
S.stat <- function(x, y, n.x=length(x), n.y=length(y), A = AUC(x, y, n.x, n.y), LT=FALSE) {
  ranks = rank(c(x,y))
  rank.x = ranks[1:n.x] ; rank.y = ranks[-1:-n.x]   # ranks within both x and y
  rank.i = rank(rank.x) ; rank.j = rank(rank.y)     # ranks within each
  S2.10 = ( sum((rank.x - rank.i)^2) - n.x*(mean(rank.x) - (n.x+1)/2)^2 ) / ((n.x-1)*n.y^2)
  S2.01 = ( sum((rank.y - rank.j)^2) - n.y*(mean(rank.y) - (n.y+1)/2)^2 ) / ((n.y-1)*n.x^2)
  if (LT) {div.LT = (A*(1-A))^2} else {div.LT=1}
  return(S = sqrt((n.x * S2.01 + n.y * S2.10)/(n.x + n.y)/div.LT))
}
V.DeLong <- function(x, y, n.x=length(x), n.y=length(y), A = AUC(x, y, n.x, n.y), LT=FALSE) {
  D2 <- D1 <- 0
  for (i in 1:n.x) {D1 = D1 + ( sum(y > x[i]) / n.y - A )^2}; D1 = D1 /(n.x * (n.x - 1))
  for (j in 1:n.y) {D2 = D2 + ( sum(y[j] > x) / n.x - A )^2}; D2 = D2 /(n.y * (n.y - 1))
  if (LT) {div.LT = (A*(1-A))^2} else {div.LT=1}
  return(V = (D1 + D2)/div.LT)
}