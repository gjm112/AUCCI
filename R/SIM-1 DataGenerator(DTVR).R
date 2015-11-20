### Part II. Simulation
## 2-1 [SIM] DataGeneration (DTVR) 2.1.1 datagenerator

## Library ############################################################################## 
library(MASS)  # for mvrnorm

## 2.1.1 datagenerator ################################################################## 
# should segregate D=0 and D=1 for R
datagenerator <- function(n, alpha0, alpha1, beta0, beta1, beta2, beta3, sig, q1, q2, q3, q4, gamma ,mu.V,Sigma, dist="Gaussian", option="VDTR"){
  V = mvrnorm(n, mu.V, Sigma)
  phi = (1+ exp(alpha0 + V%*%alpha1)^-1)^-1
  D = rbinom(n, 1, phi)
  if (option == "VD") {
    return(data.frame(disease=D, V=V))
  } else {
    mu.T = beta0 + D*beta1 + V%*%beta2 + D*V%*%beta3
    if (dist=="Gaussian") {
      Te = mu.T + sig*rnorm(n)
    } else if (dist=="Exponential"){
      Te = rexp(n,1/exp(mu.T))
    }
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