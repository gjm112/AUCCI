library(MASS)


set.seed(10)
# set.seed(10)

# should segregate D=0 and D=1 for R
datagenerator <- function(n, alpha0, alpha1, beta0, beta1, beta2, beta3, sig, q1, q2, q3, q4, gamma ,mu.V,Sigma, DT.only=FALSE){
  V = mvrnorm(n, mu.V, Sigma)
  phi = (1+ exp(alpha0 + V%*%alpha1)^-1)^-1
  D = rbinom(n, 1, phi)
  mu.T = beta0 + D*beta1 + V%*%beta2 + D*V%*%beta3
  T = mu.T + sig*rnorm(n)
  if (DT.only == FALSE) {
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
  else {
    R = 1
    DR = D
    return(data.frame(disease=D, marker=T))
  }
}

mu.V = rep(0,5)
Sigma = matrix(c(1, 0, 0.3, 0.4, -0.4, 0, 1, 0.2, 0.2, 0, 0.3, 0.2, 1, 0.7, -0.5, 0.4, 0.2, .7, 1, -0.2, -0.4, 0, -0.5, -0.2, 1), 5,5)
phi = 0.5; alpha0 = log(phi/(1-phi))
alpha1 = rep(1,5)
beta0 = 0
beta1 = 0.3484
beta2 = rep(0.1,5)
beta3 = rep(0.05,5)
sig = 1
q1 = 0.8
q2 = 0.95
q3 = 0.7
q4 = 0.90
gamma = 0.2


### theta_(1+2+..+10) almost equals average(theta_1, ... theta_10) 
{ temp <- datagenerator(10000,alpha0,alpha1,beta0,beta1,beta2,beta3,sig,gamma0,gamma1,gamma2,gamma3,mu.V,Sigma)  
  AUC(temp$marker[temp$disease==0], temp$marker[temp$disease==1])
  a <- c(1:10)
  for (i in 1:10){
    temp.temp <- temp[(10000*(i-1)+1):(10000*i),]
    a[i] <- AUC(temp.temp$marker[temp.temp$disease==0], temp.temp$marker[temp.temp$disease==1])
  }
  a
  mean(a)
}

a <- Sys.time(); temp <- datagenerator(n=100000, alpha0=alpha0, alpha1=alpha1, beta0=beta0, beta1=beta1, beta2=beta2, beta3=beta3, sig=sig, q1=q1, q2=q2, q3=q3, q4=q4, gamma=gamma, mu.V=mu.V, Sigma=Sigma, missing=FALSE) ; Sys.time()-a
a <- Sys.time(); AUC(temp$marker[temp$disease==0], temp$marker[temp$disease==1]); Sys.time()-a

# theta finder
a <- Sys.time(); n=50000000; n1=1000000; n2 = floor(n/n1)
temp <- data.frame(disease=rep(NA,n), marker=NA)
for (i in 1:n2) {temp[((i-1)*n1+1):(i*n1),] <- datagenerator(n=n1, alpha0=alpha0, alpha1=alpha1, beta0=beta0, beta1=beta1, beta2=beta2, beta3=beta3, sig=sig, q1=q1, q2=q2, q3=q3, q4=q4, gamma=gamma, mu.V=mu.V, Sigma=Sigma, DT.only=TRUE)} ; Sys.time()-a
a <- Sys.time(); AUC(temp$marker[temp$disease==0], temp$marker[temp$disease==1]); Sys.time()-a



head(temp)
Sys.time(); AUC(temp$marker[temp$disease==0], temp$marker[temp$disease==1])
Sys.time(); performance(prediction(temp$marker, temp$disease),"auc")@y.values[[1]]
Sys.time()
#temp[,c("disease","diseaseR")]
table(data.frame(disease=temp[,c("disease")], diseaseR=ifelse(is.na(temp$diseaseR),"NA",temp$disease)))

