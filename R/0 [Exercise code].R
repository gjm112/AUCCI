####### Great Greg #######
####### Pilot test codes ########
### 2-1 basic param & data  #####################################################
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
  q1 = 0.7
  q2 = 0.8
  q3 = q1
  q4 = q2
  gamma = 0.7
}
temp <- datagenerator(100,alpha0=alpha0, alpha1=alpha1, beta0=beta0, beta1=2, beta2=beta2, beta3=beta3, sig=sig, q1=q1, q2=q2, q3=q3, q4=q4, gamma=gamma, mu.V=mu.V, Sigma=Sigma, option="VDTR")


### 1-1 AUC, data2xy #####################################################
## 1.1.1
# AUC Near 1
x <- data2xy(temp)$x
y <- data2xy(temp)$y
Mee.stat(x, y, n.x=length(x), n.y=length(y))


AUC(temp$marker[temp$disease==0], temp$marker[temp$disease==1])

### 1-3 SeSp, AUC.verif #####################################################
## 1.3.1 & 1.3.2
# Example: each direct CI.method
a <- SeSp(temp,,CI.method="BG"); plot(1-a$Sp, a$Se); SeSp2AUC(a)
a <- SeSp(temp,,CI.method="MS"); plot(1-a$Sp, a$Se); SeSp2AUC(a)
a <- SeSp(temp,,CI.method="IPW"); plot(1-a$Sp, a$Se); SeSp2AUC(a)
a <- SeSp(temp,,CI.method="SP"); plot(1-a$Sp, a$Se); SeSp2AUC(a)   # some perturbations happen.
a <- SeSp(temp,,CI.method="full"); plot(1-a$Sp, a$Se); SeSp2AUC(a) 

## 1.3.3
# Example: true value #SP is closest 
AUC.verif(temp,,CI.method="BG")
AUC.verif(temp,,CI.method="MS")
AUC.verif(temp,,CI.method="IPW")
AUC.verif(temp,,CI.method="SP")
AUC.verif(temp,,CI.method="full")
AUC.verif(temp,,CI.method="naive")
AUC.verif(temp,,CI.method="He")    # basically He is equivalent to IPW

### 1-4 AUCCI.boot #####################################################
AUCCI.boot(temp,R=100,alpha=.05,base.fun=AUC.verif,CI.method="BG")
AUCCI.boot(data[1:30,], R=30, CI.method=Dir.methods, type="landscape1")

### 1-5 AUCCI.MI #####################################################
## preset arguments
m = 10  # number of imputations
MI.function = mice
MI.method = "pmm"
alpha = .05

# Example
AUCCI.MI(temp, MI.function=mice, MI.method="pmm", alpha=.05, CI.method="HM1", m = m, score.MI="fixed.r")
CI.i(temp, fun=AUCCI.MI, CI.method=CI.methods, type="landscape2",MI.function=mice, MI.method="pmm", alpha=.05, m = m)

# test
temp <- datagenerator(100,alpha0=alpha0, alpha1=alpha1, beta0=beta0, beta1=2, beta2=beta2, beta3=beta3, sig=sig, q1=q1, q2=q2, q3=q3, q4=q4, gamma=gamma, mu.V=mu.V, Sigma=Sigma, option="VDTR")
meth="NS1";
AUCCI(temp,meth,sample.var=T);AUCCI(temp,meth,LT=T)
set.seed(1);AUCCI.MI(temp,mice,"pmm", m=10, meth, sample.var=T)
set.seed(1);AUCCI.MI(temp,mice,"pmm", m=10, meth, LT=T)

## MI
complete(mice(temp,"logreg",m=5,predictorMatrix=cbind(0,(1 - diag(1, ncol(data)))[,-1]),printFlag=FALSE),2)$diseaseR
MI.norm(data)[,"diseaseR"]


### 2-1 Datagenerator #####################################################
## theta_(1+2+..+10) almost equals average(theta_1, ... theta_10) 
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
#temp[,c("disease","diseaseR")]
table(data.frame(disease=temp[,c("disease")], diseaseR=ifelse(is.na(temp$diseaseR),"NA",temp$disease)))

### 2-2 ddd #####################################################

### CI-1b optim2, multiroot.....
optim(.7, HMS2.equation, AUC.hat = 0.88, n.x=4, n.y=46, alpha=0.05, LT=FALSE, NC=FALSE, lower=0, upper=1, method="L-BFGS-B")
optim(1, HMS2.equation, AUC.hat = 0.88, n.x=4, n.y=46, alpha=0.05, LT=FALSE, NC=FALSE, lower=0, upper=1, method="L-BFGS-B")
optim2(HMS.equation, c(0.7,1), AUC.hat=0.88,n.x=4, n.y=46, alpha=0.05, LT=FALSE, NC=FALSE)
AUCCI.MI(temp.comp.error[[1]], CI.method="NS1")

multiroot(V.HM,c(.7,1.0), n.x=4, n.y=46, alpha=0.05, LT=FALSE, NC=FALSE)$root
HMS2.equation = function(...) {HMS.equation(...)^2}



### 3 Sim ########

mean(sim.data[[6]][[1]][[2]]$est.com$AUC.hat)

## Optional: saving the datafile  ########################################################
saveRDS(sim.data.212, "simdata_212.rds")
saveRDS(sim.data.311, "simdata_311.rds")
sim.data.0927 <- readRDS("sim_data_0927.rds")
saveRDS(temp.eval,"simdata_1111.rds")
## 2.3.2.2 Simulation by part (Evaluation only)  #########################################
temp.d4 = temp.d1[[1]][[1]][[1]]
a <- list
{
  # d1 = 6; d2 = 2;na.rm=ra.rm; set.seed=...
  # need to be updated if want to use this
  set.seed(100)
  for (i in 1:d1) {
    for (j in 1:d2){ 
      for (h in 1:d3) {
        temp.eval <- list()
        temp.eval[["1. complete"]] <- CI.evaluator(temp.d4[["est.com"]], param = temp.d4[["parm"]], CI.method = CI.methods, na.rm = na.rm, round=4)
        temp.eval[["2. pmm"]] <- CI.evaluator(temp.d4[["est.MI"]][["pmm"]], param = temp.d4[["parm"]], CI.method = CI.methods, na.rm = na.rm, round=4)
        temp.eval[["3. logreg"]] <- CI.evaluator(temp.d4[["est.MI"]][["logreg"]], param = temp.d4[["parm"]], CI.method = CI.methods, na.rm = na.rm, round=4)
        temp.eval[["4. imp.norm"]] <- CI.evaluator(temp.d4[["est.MI"]][["imp.norm"]], param = temp.d4[["parm"]], CI.method = CI.methods, na.rm = na.rm, round=4)
        temp.eval[["5. Bootstrap"]] <- CI.evaluator(temp.d4[["est.Dir"]], param = temp.d4[["parm"]], CI.method = Dir.methods, na.rm = na.rm, round=4)
      }
    }
  }
}


## Example: checking missing rate #####################################################
for(i in 1:20) {print(mean(sim.data[[6]][[2]][[1]][[i]]$R, na.rm=T))}
sim.data.complete <- sim.data
## addressing example
# phi = .7, theta=.8  (=> 4th)
# rho = .5            (=> 1th)
# eval                (=> 4th)    # or 5th sample data (=> [[1]][[5]])
sim.data[[4]][[1]][[4]]
head(sim.data[[4]][[1]][[1]][[5]])
head(sim.data[[4]][[1]][[3]])

t(sim.data[[4]][[2]][[4]][[1]])[,1]

