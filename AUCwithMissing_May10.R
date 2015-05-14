### log ########################################################################################
## May 6
# Q.stat : changed Q1 and Q2
## May 7 ~ 10
# truncnorm: added to library(0.1), added truncation criteria(1.1) and replace rnorm with truncnorm(2.1)
# Wilson-Score CI: added to function(0.2.4), added to sampling loop statement(2.1)
# CI summary function: extract repeating part into functions(0.3)
# Coverage probability function: extract repeating part into functions(0.4)
# added exponential samples(2.1) => but something strange...

## to do's : 
#1. check the algorithm of exponential samples(2.1)
#2. Other evaluation tools(under 0.4) : distance(th, th.hat), CI length, LNCP/RNCP, ZWI (as in Chapter 5 of draft paper)
#3. Other CI methods: NC1, NC2, CM, RG (as in Chapter 3.4 ~ 3.7)

## 0. library  ########################################################################################
# 0.1 packages
# ROCR / dplyr
library(ROCR)
library(dplyr)
library(truncnorm)  # for 2.1 sampling
library(rootSolve)  # for Newton-Raphson Method (0.2.4)

# 0.2 functions
# 0.2.1 rbivariate
rbivariate <- function(mu.x, sd.x, mu.y, sd.y, r=0.5, iter=100) {
  z1 <- rnorm(iter)
  z2 <- rnorm(iter)
  x <- sqrt(1-r^2)*sd.x*z1 + r*sd.x*z2 + mu.x
  y <- sd.y*z2 + mu.y
  return(list(x,y))
}
# 0.2.2 truncated exponential / qtexp: qunatile, rtexp: random sample of t.exp
# http://tolstoy.newcastle.edu.au/R/e10/help/10/06/8669.html
qtexp <- function(x, beta, trunc) { -log(1-x*(1-exp(-trunc/beta)))*beta }
rtexp <- function(n, beta, trunc) { qtexp(runif(n), beta, trunc) }



# 0.2.3 Q.stat : Q1 and Q2
#For Hanley-McNeil-Wald
# log: Error corrected on May 6: Q1 & Q2 switched
Q.stat <- function(data, n.x, n.y, disease="disease", marker="marker") {
  x = "$"(data, marker)["$"(data, disease)==0]
  y = "$"(data, marker)["$"(data, disease)==1]
  Q2 <- Q1 <- 0
  for (i in 1:n.x) {
    Q1 = Q1 + ( sum(y > x[i]) + sum(y==x[i])/2 )^2
  }
  for (j in 1:n.y) {
    Q2 = Q2 + ( sum(y[j] > x) + sum(y[j]==x)/2 )^2
  }
  Q1 = Q1 /(n.x * n.y^2); Q2 = Q2 /(n.x^2 * n.y)
  return(data.frame(Q1 = Q1, Q2 = Q2))
}

# 0.2.4 Wilson-Score CI
# V for variance, WS.equation for equation function, WS for CI bounds solution
V <- function(AUC, n.x, n.y, Q1 = AUC/(2-AUC), Q2=2* AUC^2 /(AUC + 1)) {
  return((AUC * (1-AUC) + (n.y -1)*(Q1 - AUC^2) + (n.x-1)*(Q2 - AUC^2))/(n.x*n.y) )
}
WS.equation <- function(AUC, theta.hat, n.x, n.y, alpha) {
  AUC + c(+1,-1)*qnorm(1-alpha/2) * sqrt(V(AUC, n.x, n.y)) - theta.hat
}
WS <- function(theta.hat, alpha = 0.05, n.x, n.y) {
  start.1 = max(0.01, theta.hat - .1); start.2 = min(0.99, theta.hat + .1)
  multiroot(WS.equation, c(start.1, start.2), theta.hat=theta.hat, n.x=n.x, n.y=n.y, alpha=alpha)$root 
}


# 0.3 CI summary function  #################################################################################
AUC.CI <- function(data, n.x, n.y, disease="disease", marker="marker", alpha) {
  z.a2 = qnorm(1-alpha/2)
  AUC.hat = (prediction(data[,marker], data[,disease]) %>% performance("auc"))@y.values[[1]]  # package ROCR
  Q1 = Q.stat(data,n.x,n.y)$Q1 ; Q2 = Q.stat(data,n.x,n.y)$Q2                   # function Q.stat HMW 
  V.HMW = V(AUC.hat, n.x, n.y, Q1, Q2)                                          # function V
  V.HM = V(AUC.hat, n.x, n.y)           # Q1, Q2 as default
  CI.WS <- matrix(NA,2,3); for (z in 1:3) (CI.WS[,z] = WS(AUC.hat, alpha[z], n.x, n.y))
  return(data.frame(AUC.hat=AUC.hat, n.x = n.x, n.y = n.y, Q1 = Q1 , Q2 = Q2 ,
                    HMW.lb = AUC.hat - z.a2*sqrt(V.HMW),  HMW.ub = AUC.hat + z.a2*sqrt(V.HMW),
                    HM.lb = AUC.hat - z.a2*sqrt(V.HM), HM.ub = AUC.hat + z.a2*sqrt(V.HM),
                    WS.lb = CI.WS[1,], WS.ub = CI.WS[2,]))
}

# 0.4 evaluation function  #################################################################################
# 0.4.1 Coverage probability 
coverage <- function(data, dim = c(5,9,3), CI.type=c("HMW", "HM", "WS")) {
  # data: list of inferences(AUC hat, CI's lower/upper bounds)
  coverage <- array(,c(dim,length(CI.type)))
  for (z in 1:length(CI.type)) {
    CI.lb = paste0(CI.type[z],".lb") ; CI.ub = paste0(CI.type[z],".ub")
    for (i in 1:dim[1]) {        # AUC (0.6 ~ 0.95)
      for (j in 1:dim[2]){       # n*pi (30,50,200) * (0.3, 0.5, 0.7)
        for (k in 1:dim[3]){     # alpha 10%, 5%, 1%
          coverage[i,j,k, z] <- mean(sapply(data[[i]][[j]],function(x) {
            # x$HMW.lb[k] <= x$AUC & x$HMW.ub[k] >= x$AUC }))       # replace this!!!
            
            x[k,CI.lb] <= x$AUC & x[k, CI.ub] >= x$AUC }), na.rm=T)       # with this..
        }
      }
    }
  }
  cat(" row: AUC = 0.6 ~ 0.95 \n col: n*pi = (30,50,200)*(.3,.5,.7) \n 3rd: alpha = .1, .05, .01 \n 4th: CI methods = HMW, HM, WS")
  return(coverage)
}



## 1. population data  ######################################################################################################
## 1.1 ~1.2 parameters
## 1.1 parms1 (theta, mu, sig, beta)        : dim(parm) = 5AUC + z.a2*sqrt(V.HM),

# theta = AUC
parm1 <- data.frame(theta = c(0.6, 0.7, 0.8, 0.9, 0.95))
parm1 <- within(parm1,
{ # Under normal distribution
  mu.x <- 1; sig.x <- 1; sig.y <- 2
  # Under exponential distribution
  beta.x <- 1
  # Residual params under thetas
  mu.y <- mu.x + qnorm(theta)*sqrt(sig.x^2 + sig.y^2)
  beta.y <- theta * beta.x / (1-theta)
  # bounds for trunction
  a.x <- mu.x - 3* sig.x; b.x <- mu.x + 3*sig.x   # bounds for truncated noraml
  a.y <- mu.y - 3* sig.y; b.y <- mu.y + 3*sig.y
  c.x <- qexp(0.99, beta.x); c.y <- qexp(0.99, beta.y) # bounds for truncated exp
})

## 1.2 parms2: n, pi
# n= sample size,  pi = prevalence rate     : dim(parm2) = 3*3=9
# n.x and n.y (= n * pi) are assumed to be fixed not random.
parm2 <- data.frame(n = rep(c(30, 50, 200), each=3), pi = rep(c(0.30, 0.50, 0.70),3))


## 1.3 alpha levels
alpha <- c(0.1, 0.05, 0.01)
z.a2<- qnorm(1-alpha/2)     # z_alpha/2



## 2. Sample data  ##################################################################################################
## 2.1 generating sample data with bivariate-normal distribution

sample.normal <- sample.normal.stat <- list()   # outer shell   (data & statistics) normal
temp.2 <- temp.2.stat <- list()                 # middle shells
temp.1 <- temp.1.stat <- list()                 # inner shells
sample.exp <- sample.exp.stat <- list()         # outer shell   (data & statistics) exponential
temp.2.exp <- temp.2.exp.stat <- list()         # middle shells
temp.1.exp <- temp.1.exp.stat <- list()         # inner shells
{
  set.seed(10)
  for (i in 1:5) { # i: AUC
    for (j in 1:9) { # j: n sample size, pi prevalence rate
      mu.x <- parm1$mu.x[i]; mu.y <- parm1$mu.y[i]; sig.x <- parm1$sig.x[i]; sig.y <- parm1$sig.y[i];
      beta.x <- parm1$beta.x[i]; beta.y <- parm1$beta.y[i]; theta <- parm1$theta[i]
      a.x <- parm1$a.x[i]; b.x <- parm1$b.x[i]; a.y <- parm1$a.y[i]; b.y <- parm1$b.y[i]
      c.x <- parm1$c.x[i]; c.y <- parm1$c.y[i]
      n <- parm2$n[j] ; pi <- parm2$pi[j]
      n.y <- n*pi; n.x <- n - n.y
      for (k in 1:10000){ # k: num of simulations(samples)
        # normal
        temp.1[[k]] <- temp <- data.frame(disease = c(rep(0,n.x), rep(1, n.y)),
                                          marker  = c(rtruncnorm(n-n*pi,a.x, b.x, mu.x, sig.x), rtruncnorm(n*pi,a.y, b.y, mu.y, sig.y)))
        temp.1.stat[[k]] <- AUC.CI(temp, n.x, n.y, alpha=alpha)
        temp.1.stat[[k]]$AUC <- theta
        # exponential
        temp.1.exp[[k]] <- temp <- data.frame(disease = c(rep(0,n.x), rep(1, n.y)),
                                              marker  = c(rtexp(n-n*pi, beta.x, c.x), rtexp(n*pi,beta.y, c.y)))
        temp.1.exp.stat[[k]] <- AUC.CI(temp, n.x, n.y, alpha=alpha)
        temp.1.exp.stat[[k]]$AUC <- theta
      }
      temp.2[[j]] <- temp.1
      temp.2.stat[[j]] <- temp.1.stat
      temp.2.exp[[j]] <- temp.1.exp
      temp.2.exp.stat[[j]] <- temp.1.exp.stat
    }
    sample.normal[[i]] <- temp.2
    sample.normal.stat[[i]] <- temp.2.stat
    sample.exp[[i]] <- temp.2.exp
    sample.exp.stat[[i]] <- temp.2.exp.stat
  }
}


## 3. Evaluation  ########################################################################################
#3.1 Coverage probability

print(coverage.normal <- coverage(sample.normal.stat))
print(coverage.exp <- coverage(sample.exp.stat))
