### Part II. Simulation
## 2-1 [SIM] DataGeneration (DTVR) 2.1.1 datagenerator

## Library ############################################################################## 
library(MASS)  # for mvrnorm

## 2.1.1 datagenerator ################################################################## 
# should segregate D=0 and D=1 for R
datagenerator <- function(n, alpha0, alpha1, beta0, beta1, beta2, beta3, sig, q1, q2, q3, q4, gamma ,mu.V,Sigma, option="VDTR"){
  V = mvrnorm(n, mu.V, Sigma)
  phi = (1+ exp(alpha0 + V%*%alpha1)^-1)^-1
  D = rbinom(n, 1, phi)
  if (option == "VD") {
    return(data.frame(disease=D, V=V))
  } else {
    mu.T = beta0 + D*beta1 + V%*%beta2 + D*V%*%beta3
    Te = mu.T + sig*rnorm(n)

    if (option == "VDR") {
      return(data.frame(disease=D, marker=Te))
    } else {  # option="VDTR"
      # for R(missingness)
      R.T = (Te <= quantile(Te, q1))*(1-D) + (Te <= quantile(Te, q3))*D
      Q2 = apply(V, 2, quantile, q2)
      Q4 = apply(V, 2, quantile, q4)
      normalV = t(apply(V,1,"<=",Q2))*(1-D) + t(apply(V,1,"<=",Q4))*D
      R.V = as.vector(t(apply(normalV,1,prod)))
      rho = (R.T*R.V)*gamma
      R = rbinom(n, 1, rho)
      DR = ifelse(R==0, D, NA)    
      return(data.frame(disease=D, marker=Te, V=V, R=R, diseaseR=DR))
    }
  }
}

## parameters for pilot test 
{
  mu.V = rep(0,5)
  Sigma = matrix(c(1, 0, 0.3, 0.4, -0.4, 0, 1, 0.2, 0.2, 0, 0.3, 0.2, 1, 0.7, -0.5, 0.4, 0.2, .7, 1, -0.2, -0.4, 0, -0.5, -0.2, 1), 5,5)
  alpha0 = 0
  alpha1 = rep(1,5)
  beta0 = 0
  beta1 = 0.8089
  beta2 = rep(0.1,5)
  beta3 = rep(0.05,5)
  sig = 1
  q1 = 0.8
  q2 = 0.95
  q3 = 0.7
  q4 = 0.90
  gamma = 0.2
}

### theta_(1+2+..+10) almost equals average(theta_1, ... theta_10) 
{ temp <- datagenerator(10000, alpha0=alpha0, alpha1=alpha1, beta0=beta0, beta1=beta1, beta2=beta2, beta3=beta3, sig=sig, q1=q1, q2=q2, q3=q3, q4=q4, gamma=gamma, mu.V=mu.V, Sigma=Sigma, option="VDTR")  
  AUC(temp$marker[temp$disease==0], temp$marker[temp$disease==1])
  a <- c(1:10)
  for (i in 1:10){
    temp.temp <- temp[(10000*(i-1)+1):(10000*i),]
    a[i] <- AUC(temp.temp$marker[temp.temp$disease==0], temp.temp$marker[temp.temp$disease==1])
  }
  a
  mean(a)
}


Sys.time(); AUC(temp$marker[temp$disease==0], temp$marker[temp$disease==1])

#temp[,c("disease","diseaseR")]
table(data.frame(disease=temp[,c("disease")], diseaseR=ifelse(is.na(temp$diseaseR),"NA",temp$disease)))