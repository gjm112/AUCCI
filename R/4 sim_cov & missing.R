library(MASS)


set.seed(10)
# set.seed(10)


datagenerator <- function(n,alpha0,alpha1,beta0,beta1,beta2,beta3,sig,gamma0,gamma1,gamma2,gamma3,mu.V,Sigma){
  V = mvrnorm(n, mu.V, Sigma)
  phi = (1+ exp(alpha0 + V%*%alpha1)^-1)^-1
  D = rbinom(n, 1, phi)
  mu.T = beta0 + D*beta1 + V%*%beta2 + D*V%*%beta3
  T = mu.T + sig*rnorm(n)
  rho = (1+ exp(gamma0 + D*gamma1 + T*gamma2 + V%*%gamma3)^-1)^-1
  R = rbinom(n, 1, rho)
  DR = ifelse(R==0, D, NA)
  return(data.frame(disease=D, marker=T, V=V, R=R, diseaseR=DR))
}

mu.V = rep(0,5)
Sigma = matrix(c(1, 0, 0.3, 0.4, -0.4, 0, 1, 0.2, 0.2, 0, 0.3, 0.2, 1, 0.7, -0.5, 0.4, 0.2, .7, 1, -0.2, -0.4, 0, -0.5, -0.2, 1), 5,5)
phi = 0.3; alpha0 = log(phi/(1-phi))
alpha1 = rep(1,5)
beta0 = 0
beta1 = 0.1
beta2 = rep(0.1,5)
beta3 = rep(0.05,5)
sig = 1
gamma0 = -1
gamma1 = -0.2
gamma2 = rep(-0.1,5)
gamma3 = rep(-0.05,5)

set.seed(1)
temp <- datagenerator(100000,alpha0,alpha1,beta0,beta1,beta2,beta3,sig,gamma0,gamma1,gamma2,gamma3,mu.V,Sigma)
AUC(temp$marker[temp$disease==0], temp$marker[temp$disease==1])

#temp[,c("disease","diseaseR")]
#table(data.frame(disease=temp[,c("disease")], diseaseR=ifelse(is.na(temp$diseaseR),"NA",temp$disease)))

