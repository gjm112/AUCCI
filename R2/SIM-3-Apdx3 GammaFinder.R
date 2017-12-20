source("R/CI-1 Base functions.R")
source("R/CI-1b Auxiliary functions.R")
source("R/CI-2 AUCCI (Complete).R")
source("R/SIM-1 DataGenerator(DTVR).R")

# parameters (except alpha0 and beta1)
mu.V = rep(0,5)
Sigma = matrix(c(1, 0, 0.3, 0.4, -0.4, 0, 1, 0.2, 0.2, 0, 0.3, 0.2, 1, 0.7, -0.5, 0.4, 0.2, .7, 1, -0.2, -0.4, 0, -0.5, -0.2, 1), 5,5)
alpha1 = rep(1,5)
beta0 = 0
beta2 = rep(0.1,5)
beta3 = rep(0.05,5)
sig = 1

alpha0 = 0
beta1 = c(0.8089, 1.4486, 1.97674, 2.96704, .8319, 1.47286, 2.00192, 2.9939)
q1 <- q3 <- .7; q2 <- q4 <- .8; gamma <- .7

meanR0 <- function(alpha0, beta1, q1, q2, gamma) {
  q3 <- q1; q4 <- q2
  return(mean(datagenerator(n=1000, alpha0=alpha0, alpha1=alpha1, beta0=beta0, beta1=beta1, beta2=beta2, beta3=beta3, sig=sig, q1=q1, q2=q2, q3=q3, q4=q4, gamma=gamma, mu.V=mu.V, Sigma=Sigma, dist="Gaussian", option="VDTR")$R) )
}
meanR1 <- Vectorize(meanR0, c("alpha0", "beta1"))
meanR2 <- function(q1, q2, gamma) {
  result <- meanR1(alpha0 = rep(c(0,1.6111), each=4), beta1 = beta1, q1 = q1, q2 = q2, gamma =gamma)  
  result <- c(mean=mean(result), result)
  return(result)
}
meanR2(.7, .8, .7)

# meanR3 for grid search
meanR3 <- Vectorize(meanR2, c("q1", "q2", "gamma"))
grid.tmp <- expand.grid(q1 = (1:5)*.05 + .7, q2 = (1:5)*.03 + .8, gamma = (1:5)*.05 + .7)
set.seed(10)
a <- meanR3(grid.tmp[,1], grid.tmp[,2], grid.tmp[,3])
a[1,]

set.seed(10)
meanR2(.85, .9, .9)  #.51
meanR2(.9, .95, .95)  #.70
meanR2(.99, .99, .95)  #.90

meanR2(.7, .8, .7)  #.22
meanR2(.8, .9, .8)  #.43
