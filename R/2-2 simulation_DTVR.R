library(MASS)  # for mvrnorm


set.seed(10)

# should segregate D=0 and D=1 for R
datagenerator <- function(n, alpha0, alpha1, beta0, beta1, beta2, beta3, sig, q1, q2, q3, q4, gamma ,mu.V,Sigma, option="VDTR"){
  V = mvrnorm(n, mu.V, Sigma)
  phi = (1+ exp(alpha0 + V%*%alpha1)^-1)^-1
  D = rbinom(n, 1, phi)
  if (option == "VD") {
    return(data.frame(disease=D, V=V))
  } else {
    mu.T = beta0 + D*beta1 + V%*%beta2 + D*V%*%beta3
    T = mu.T + sig*rnorm(n)

    if (option == "VDR") {
      return(data.frame(disease=D, marker=T))
    } else {  # option="VDTR"
      # for R(missingness)
      R.T = (T <= quantile(T, q1))*(1-D) + (T <= quantile(T, q3))*D
      Q2 = apply(V, 2, quantile, q2)
      Q4 = apply(V, 2, quantile, q4)
      normalV = t(apply(V,1,"<=",Q2))*(1-D) + t(apply(V,1,"<=",Q4))*D
      R.V = as.vector(t(apply(normalV,1,prod)))
      rho = (R.T*R.V)*gamma
      R = rbinom(n, 1, rho)
      DR = ifelse(R==0, D, NA)    
      return(data.frame(disease=D, marker=T, V=V, R=R, diseaseR=DR))
    }
  }
}
parm.set <- list()

mu.V = rep(0,5)
Sigma = matrix(c(1, 0, 0.3, 0.4, -0.4, 0, 1, 0.2, 0.2, 0, 0.3, 0.2, 1, 0.7, -0.5, 0.4, 0.2, .7, 1, -0.2, -0.4, 0, -0.5, -0.2, 1), 5,5)
alpha0 = 1.6111   # for phi=0.7
alpha1 = rep(1,5)
beta0 = 0
### beta1 = 0.9 below->
beta2 = rep(0.1,5)
beta3 = rep(0.05,5)
sig = 1
q1 = 0.8
q2 = 0.95
q3 = 0.7
q4 = 0.90
gamma = 0.2


### theta_(1+2+..+10) almost equals average(theta_1, ... theta_10) 
{ temp <- datagenerator(10000, alpha0=alpha0, alpha1=alpha1, beta0=beta0, beta1=beta1, beta2=beta2, beta3=beta3, sig=sig, q1=q1, q2=q2, q3=q3, q4=q4, gamma=gamma, mu.V=mu.V, Sigma=Sigma, option="VDT")  
  AUC(temp$marker[temp$disease==0], temp$marker[temp$disease==1])
  a <- c(1:10)
  for (i in 1:10){
    temp.temp <- temp[(10000*(i-1)+1):(10000*i),]
    a[i] <- AUC(temp.temp$marker[temp.temp$disease==0], temp.temp$marker[temp.temp$disease==1])
  }
  a
  mean(a)
}

### free parameters: alpha0(for phi), beta1(for theta)
### To estimate parameters that satisfy specific phi and theta,
### refer to 2-2-1 parameterfinder.R(beta1) and 2-2-2 alphafinder.R(alpha0).
alphabet = data.frame(phi = rep(c(.5,.7), each=3), theta = rep(c(.7,.8,.9),2), alpha0 = rep(c(0,1.6111), each=3), beta1=rep(NA,6))
alphabet$beta1 = c(0.8089, 1.4486, 1.97674, NA, NA, NA)


head(temp)
Sys.time(); AUC(temp$marker[temp$disease==0], temp$marker[temp$disease==1])

#temp[,c("disease","diseaseR")]
table(data.frame(disease=temp[,c("disease")], diseaseR=ifelse(is.na(temp$diseaseR),"NA",temp$disease)))

