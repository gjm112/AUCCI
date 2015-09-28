### Part I. CI functions
## 1-1 [CI] Base functions: 1.1.1 CI.base, 1.1.2 Variance estimators

# Library ##############################################################################
library(rootSolve)  # for Newton-Raphson Method (0.2.4)
library(sn)         # for Owen's T-function in Double Gaussian Method (0.2.15)

# CI methods :    ref to 1.CI test(temp), 2.simulation, 3.evaluation
CI.methods = c("HanleyMcNeilWald", "HanleyMcNeilExponential", "HanleyMcNeilScore", 
               "NewcombeExponential", "NewcombeScore", "CortesMohri", "ReiserGuttman", 
               "ClopperPearson", "Bamber", "WilsonScore",
               "HalperinMee", "DoubleBeta", "DoubleGaussian",
               "MannWhitney", "MannWhitneyLT", "DeLong")
########################################################################################

## 1.1.1 Basic internal functions ######################################################
# CI.base, AUC.classic, AUC
CI.base <- function(mu.hat, Var, alpha) {
  return(mu.hat + c(-1,1) * qnorm(1-alpha/2)*sqrt(Var))
}
AUC.classic <- function(x, y, n.x=length(x), n.y=length(y)) {
  A <- 0; for (i in 1:n.x) { A = A + ( sum(y > x[i]) + sum(y==x[i])/2 ) }
  return(A = A/ n.x / n.y)
}
# new AUC function added with much improved efficiency
AUC <- function(x, y, n.x=length(x), n.y=length(y)) {
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

## 1.1.2 variance estimators ###########################################################
# Q.stat for HM and its derivatives
Q.stat <- function(x, y, n.x=length(x), n.y=length(y)) {
  Q2 <- Q1 <- 0
  for (i in 1:n.x) {Q1 = Q1 + ( sum(y > x[i]) + sum(y==x[i])/2 )^2}; Q1 = Q1 /(n.x * n.y^2)
  for (j in 1:n.y) {Q2 = Q2 + ( sum(y[j] > x) + sum(y[j]==x)/2 )^2}; Q2 = Q2 /(n.x^2 * n.y)
  return(data.frame(Q1 = Q1, Q2 = Q2))
}
V.HM <- function(AUC, n.x, n.y, Q1 = AUC/(2-AUC), Q2=2* AUC^2 /(AUC + 1)) {
  # default: HME method
  return((AUC * (1-AUC) + (n.y -1)*(Q1 - AUC^2) + (n.x-1)*(Q2 - AUC^2))/(n.x*n.y) )
}
HMS.equation = function(AUC, AUC.hat, n.x, n.y, alpha) {
  AUC + c(+1,-1)*qnorm(1-alpha/2) * sqrt(V.HM(AUC, n.x, n.y)) - AUC.hat
}
V.NC <- function(AUC, n.x, n.y) {
  N = (n.x+n.y)/2
  return( AUC*(1-AUC)/(n.x*n.y) * ( (2*N-1) - (3*N-3) / ((2-AUC)*(AUC+1)) ) )
}
NC.equation <- function(AUC, AUC.hat, n.x, n.y, alpha) {
  AUC + c(+1,-1)*qnorm(1-alpha/2) * sqrt(V.NC(AUC, n.x, n.y)) - AUC.hat   
}
V.CM <- function(AUC, n.x, n.y){
  z1 = function (w, n.x, n.y, k) {1-0.5*(w/n.x+(k-w)/n.y)}
  z2 = function (w, n.x, n.y, k) 
    (n.y*w^2+n.x*(k-w)^2+n.y*(n.y+1)*w+n.x*(n.x+1)*(k-w)-2*w*(k-w)*(n.y+n.x+1))/(12*n.y^2*n.x^2)
  F.base <- function (w, n.x, n.y, k) {choose(n.y-k+2*w,w)*choose(n.x+k-2*w,k-w)}
  k=round((n.x+n.y)*(1-AUC))
  F1=F2=F3=F.den=0
  for (w in 0:k){
    F1 = F1 + F.base(w, n.x, n.y, k)*z1(w, n.x, n.y, k)^2
    F.den = F.den + F.base(w, n.x, n.y, k)
    F2 = F2 + F.base(w, n.x, n.y, k)*z1(w, n.x, n.y, k)
    F3 = F3 + F.base(w, n.x, n.y, k)*z2(w, n.x, n.y, k)
  }
  return(V = (F1/F.den - F2^2/F.den^2 + F3/F.den))
}
V.Bm <- function(x, y, n.x=length(x), n.y=length(y), AUC.hat=AUC(x,y,n.x,n.y)){
  b.yyx <- b.xxy <- p.xy <- 0
  for (i in 1:n.x) {b.yyx = b.yyx + sum(y < x[i])^2 + sum(x[i] < y)^2 -2*sum(y < x[i])*sum(x[i] < y) }
  b.yyx = b.yyx/(n.x*n.y^2)
  for (j in 1:n.y) {b.xxy = b.xxy + sum(x < y[j])^2 + sum(y[j] < x)^2 -2*sum(x < y[j])*sum(y[j] < x) }
  b.xxy = b.xxy/(n.y*n.x^2)
  for (i in 1:n.x) {  p.xy = p.xy + sum(y != x[i] ) }
  p.xy = p.xy /(n.x*n.y)
  return(1/(4*(n.y-1)*(n.x-1))*(p.xy + (n.y -1)*b.xxy + (n.x -1)*b.yyx - 4*(n.y + n.x -1)*(AUC.hat -0.5)^2))
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
Halperin.stat <- function(x, y, n.x=length(x), n.y=length(y)) {
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
Halperin.equation <- function(AUC, AUC.hat, N.J.hat, alpha){
  AUC + c(+1,-1)*qnorm(1-alpha/2) * sqrt(AUC*(1-AUC)/N.J.hat) - AUC.hat
}
Halperin <- function(x, y, n.x=length(x), n.y=length(y), alpha) {
  Halperin.stat = Halperin.stat(x, y, n.x, n.y)
  AUC.hat = Halperin.stat$AUC.hat
  N.J.hat = Halperin.stat$N.J.hat
  multiroot(Halperin.equation, c(0.01, 0.99), AUC.hat=AUC.hat, N.J.hat=N.J.hat, alpha=alpha)$root 
}
V.DB <- function(alp, n.x, n.y) {
  R1 = gamma(alp + 1)^2 / gamma(2*alp + 1)
  R2 = gamma(2*alp +1)*gamma(alp +1)/gamma(3*alp+1)
  theta = 1 - R1
  Q = 1 - 2*R1 + R2
  V = (theta*(1-theta) + (n.x+n.y-2)*(Q-theta^2))/n.x/n.y
}
DB.equation <- function(alp, AUC.hat, n.x, n.y, alpha) {
  AUC = 1 - gamma(alp + 1)^2 / gamma(2*alp + 1)
  return(AUC + c(+1,-1)*qnorm(1-alpha/2) * sqrt(V.DB(alp, n.x, n.y)) - AUC.hat)
}
DB <- function(AUC.hat, n.x, n.y, alpha) {
  alp <- multiroot(DB.equation, c(1,3), AUC.hat=AUC.hat, n.x=n.x, n.y=n.y, alpha=alpha)$root
  return(1 - gamma(alp + 1)^2 / gamma(2*alp + 1))
}
V.DG <- function(delta, n.x, n.y) {
  theta = pnorm(delta/sqrt(2))
  Owen = T.Owen(delta/sqrt(2), 1/sqrt(3))
  return(V = ((n.x + n.y - 1)*theta*(1-theta) - 2*(n.x + n.y -2) * Owen ) / (n.x * n.y))
}
DG.equation <- function(delta, AUC.hat, n.x, n.y, alpha) {
  AUC = pnorm(delta/sqrt(2))
  return(AUC + c(+1,-1)*qnorm(1-alpha/2) * sqrt(V.DG(delta, n.x, n.y)) - AUC.hat)
}
DG <- function(AUC.hat, n.x, n.y, alpha) {
  delta <- multiroot(DG.equation, c(1, 3), AUC.hat=AUC.hat, n.x=n.x, n.y=n.y, alpha=alpha)$root
  return(pnorm(delta/sqrt(2)))
}
S.stat <- function(x, y, n.x=length(x), n.y=length(y)) {
  ranks = rank(c(x,y))
  rank.x = ranks[1:n.x] ; rank.y = ranks[-1:-n.x]   # ranks within both x and y
  rank.i = rank(rank.x) ; rank.j = rank(rank.y)     # ranks within each
  S2.10 = ( sum((rank.x - rank.i)^2) - n.x*(mean(rank.x) - (n.x+1)/2)^2 ) / ((n.x-1)*n.y^2)
  S2.01 = ( sum((rank.y - rank.j)^2) - n.y*(mean(rank.y) - (n.y+1)/2)^2 ) / ((n.y-1)*n.x^2)
  return(S = sqrt((n.x * S2.01 + n.y * S2.10)/(n.x + n.y)))
}
V.DeLong <- function(x, y, n.x=length(x), n.y=length(y), A = AUC(x, y, n.x, n.y)) {
  D2 <- D1 <- 0
  for (i in 1:n.x) {D1 = D1 + ( sum(y > x[i]) / n.y - A )^2}; D1 = D1 /(n.x * (n.x - 1))
  for (j in 1:n.y) {D2 = D2 + ( sum(y[j] > x) / n.x - A )^2}; D2 = D2 /(n.y * (n.y - 1))
  return(V = D1 + D2)
}