### Part II. Simulation
## 2-3 [SIM] Simulator:  2.3.1 Parameters,  2.3.2 Simulator

### !!! must replace 5th, 6th elements of param beta1, theta, theta.SE!!!!


## 2.3.1. Parameters  ###################################################################

## 2.3.1.1 Dimension definition #########################################################
{
  # 1. phi x theta (6)
  # 1.1 phi        (2)         .5, .7
  # 1.2 theta      (3)         .8, .9, .95  
  # 2. rho        (2)         approx .25, .5
  # 3. n          (3)         200, 100, 50
  # 4. CI.method  (54)
  # 4.1 Combined methods(48)
  # 4.1.1 MI.methods (3)            MICE(pmm, logistic), NORM
  # 4.1.2 CI.methods (16)           HMW, HME, ...
  # 4.2 Direct CI.method (6)     Naive, BG, MS, IPW, SP, He
}

## 2.3.1.2 Param1 - Distributional #######################################################
mu.V = rep(0,5)
Sigma = matrix(c(1, 0, 0.3, 0.4, -0.4, 0, 1, 0.2, 0.2, 0, 0.3, 0.2, 1, 0.7, -0.5, 0.4, 0.2, .7, 1, -0.2, -0.4, 0, -0.5, -0.2, 1), 5,5)
# alpha0
alpha1 = rep(1,5)
beta0 = 0
# beta1
beta2 = rep(0.1,5)
beta3 = rep(0.05,5)
sig = 1

# free parameters: alpha0(for phi), beta1(for theta)
# To estimate parameters that satisfy specific phi and theta,
# refer to 2-2-1 parameterfinder.R(beta1) and 2-2-2 alphafinder.R(alpha0).
param1 = list(alphabet = data.frame(phi = rep(c(.5,.7), each=3), 
                                    theta = c(0.7999914248, 0.89999851, 0.950002331, 0.800002545, 0.899999218, 0.950000849), 
                                    theta.SE =c(5.80868E-06,4.09076E-06,2.76121E-06,7.79383E-06, 5.34159E-06, 3.42365E-06) ,
                                    alpha0 = rep(c(0,1.6111), each=3), 
                                    beta1=c(0.8089, 1.4486, 1.97674,.8319,1.47286,2.00192)), 
              gamma = data.frame(q1 = c(.7, .8), q2 = c(.8, .9), q3=c(.7, .8), q4=c(.8, .9), gamma = c(.7, .8)))
# alpha0 = 0 for phi==0, alpha0 = 1.6111 for phi==.7
# q3= q1, q4= q2: MAR settings

## 2.3.1.3 Param2 - Simulation ###########################################################
n = c(200, 100, 50); n.max = max(n)
n.sim = 10000    # To be changed to 10,000
alpha = 0.05

# CI.methods from 1-1 [CI] Base functions.R
CI.methods = c("Bm", "HM1", "HM2", "NS1", "NW", "NS2", "Mee", "DL", "RG", "DB", "DG", "CM")

# MI.methods from 1-5 [CI] AUCCI_MI
MI.methods = data.frame(functions=c(rep("mice",2),rep("norm",1)), methods=c("pmm","logreg","imp.norm"))
m = 10

# Direct methods with bootstrap. R = # of resamples
Dir.methods = c("naive", "BG", "MS", "SP", "IPW", "He")
R = 30

## 2.3.2 Simulator #######################################################################

## 2.3.2.1 Simulation package (Data generation + Inference + Evaluation)  ################
set.seed(11); n.sim=1000
{ 
  # setting dimensions and seed
  d1 = 1:length(param1$alphabet$theta); d2 = 1:dim(param1$gamma)[1]; d3 = 1:length(n); set.seed = 50
  imputation=TRUE
  # d1 <- 2; d2 <-1; d3 <-2 ; imputation=FALSE
  d123 <- length(d1)*length(d2)*length(d3)
  na.rm = TRUE # for evaluation
  
  # data structure
  # temp.d1 (=sim.data) = list of 6 elements(temp.d2): one for each theta & phi
  # temp.d2             = list of 2 elements(unnamed list): one for each rho
  # temp.d3             = list of 4 elements(temp.sim, param, temp.est, temp.eval)
  # temp.sim            = list of 10000(=n.sim) elements: the unit samples
  # temp.est            = list of 3 elements(temp.est.n): one for n(200, 100, 5)
  # temp.est.n          = a dataframe of estimators(point and CI's)
  # temp.eval           = a dataframe of evaluation of the estimators
  
  # creating empty lists
  sim.data  <- temp.d1 <- temp.d2 <- temp.d3 <- temp.sim <- temp.n <- temp <- list() 
    
  # loop
  for (i in d1) {
    alpha0 = param1$alphabet$alpha0[i]
    beta1  = param1$alphabet$beta1[i]
    phi    = param1$alphabet$phi[i]
    theta  = param1$alphabet$theta[i]
    
    for (j in d2){
      gamma  = param1$gamma$gamma[j]
      q1  = param1$gamma$q1[j]
      q2  = param1$gamma$q2[j]
      q3  = param1$gamma$q3[j]
      q4  = param1$gamma$q4[j]
        
      for (h in d3) {
        pb <- txtProgressBar(min=0, max = n.sim, char = paste0(((i-1)*2+j-1)*3+h, "/",d123," (DataGen)|"), style=3)

        # A. Data generation(from 2.)
        for (k in 1:n.sim) {
          # temp.sim[[k]] <- temp <- datagenerator.normal(n=n[h], theta=theta, phi=phi, sig.y=3)
          temp.sim[[k]] <- temp <- datagenerator(n=n[h], alpha0=alpha0, alpha1=alpha1, beta0=beta0, beta1=beta1, beta2=beta2, beta3=beta3, sig=sig, q1=q1, q2=q2, q3=q3, q4=q4, gamma=gamma, mu.V=mu.V, Sigma=Sigma, option="VDTR")
          setTxtProgressBar(pb,k)
        }
        temp.n[[h]] <- temp.sim
        
        # B. Imputation
        if (imputation==TRUE) {
        # input: temp.n[[h]] / output: temp.MI
        temp.MI <- temp.MI.sim <- list()
        for (l in MI.methods[,"methods"]) {
          pb <- txtProgressBar(min=0, max = n.sim, char = paste0(((i-1)*2+j-1)*3+h, "/",d123," (MI-",l,")|"), style=3)
          if (l == "pmm" |l == "logreg") {
            for (k in 1:n.sim) {
              temp.imp <- mice(data=temp.sim[[k]], method=l, m=m, predictorMatrix=cbind(0,(1 - diag(1, ncol(temp.sim[[k]])))[,-1]), printFlag=FALSE)
              temp.comp <- list()
              for (q in 1:m) {temp.comp[[q]] = complete(temp.imp, action = q)[,c("diseaseR", "marker")]}
              temp.MI.sim[[k]] <- temp.comp
              setTxtProgressBar(pb,k)
            }
            temp.MI[[l]] = temp.MI.sim
          }
          else if (l == "imp.norm") {
            ###TBD!!!!!!
            temp.MI[[l]] = temp.MI[[1]]  # proxy (pmm is copied)
          }
          else {print("Wrong MI method!")}
        }
        } #imputation==TRUE
        
        # C-1. Inferences (complete datasets)
        pb <- txtProgressBar(min=0, max = n.sim, char = paste0(((i-1)*2+j-1)*3+h, "/",d123," (Est.com-)|"), style=3)  
        temp.est <- as.data.frame(matrix(NA,n.sim,(length(CI.methods)*2+1)))
        names(temp.est) <- c("AUC.hat",paste0(rep(CI.methods,each=2),c(".lb",".ub")))
        for (k in 1:n.sim) {
          temp.est[k,] <- CI.i(temp.sim[[k]],fun=AUCCI, CI.method=CI.methods, alpha=alpha, type="landscape2")
          setTxtProgressBar(pb,k)
        }
        
        # C-2. Inferences (incomplete datasets with MI)
        if (imputation==TRUE) {
        temp.estMI <- list()
        temp.estMI.n <- as.data.frame(matrix(NA,n.sim,(length(CI.methods)*2+1)))
        names(temp.estMI.n) <- c("AUC.hat",paste0(rep(CI.methods,each=2),c(".lb",".ub")))
        for (l in MI.methods[,"methods"]) {
          pb <- txtProgressBar(min=0, max = n.sim*3, char = paste0(((i-1)*2+j-1)*3+h, "/",d123," (Est.MI-",l,")|"), style=3)
          l.index = which(l ==MI.methods[,"methods"] )
          for (k in 1:n.sim) {
            temp.estMI.n[k,] <- CI.i(temp.MI[[l]][[k]], fun=AUCCI.MI, CI.method=CI.methods, m=m, alpha=alpha, type="landscape2")
            setTxtProgressBar(pb,k + (l.index-1)*n.sim)
          }
          temp.estMI[[l]] <- temp.estMI.n
        }
        
        # C-3. Inferences (incomplete datasets with Direct methods(Bootstrapping))
        temp.estDir <- as.data.frame(matrix(NA,n.sim,(length(Dir.methods)*2+1)))
        names(temp.estDir) <- c("AUC.hat",paste0(rep(Dir.methods,each=2),c(".lb",".ub")))
        pb <- txtProgressBar(min=0, max = n.sim, char = paste0(((i-1)*2+j-1)*3+h, "/",d123," (Est.Dir-)|"), style=3)
        for (k in 1:n.sim) {
          temp.estDir[k,] <- CI.i(data=temp.n[[h]][[k]], fun=AUCCI.boot, R=R, CI.method=Dir.methods, alpha=alpha, type="landscape2")
          setTxtProgressBar(pb,k)
        }
        } #imputation==TRUE
        
        # storing A ~ C-2
        
        temp.d4 <- list(parm = data.frame(theta=theta, phi=phi, gamma=gamma),
                                        sim = temp.sim, est.com = temp.est)
        if (imputation==TRUE) {temp.d4[["MI"]] = temp.MI
                               temp.d4[["est.MI"]] = temp.estMI
                               temp.d4[["est.Dir"]] = temp.estDir} 
        
        # D-1. Evaluation of C-1(est.com)
        temp.eval <- list()
        temp.eval[["1. complete"]] <- CI.evaluator(temp.d4[["est.com"]], param = temp.d4[["parm"]], CI.method = CI.methods, na.rm = na.rm, round=4)
        if (imputation==TRUE) {
        temp.eval[["2. pmm"]] <- CI.evaluator(temp.d4[["est.MI"]][["pmm"]], param = temp.d4[["parm"]], CI.method = CI.methods, na.rm = na.rm, round=4)
        temp.eval[["3. logreg"]] <- CI.evaluator(temp.d4[["est.MI"]][["logreg"]], param = temp.d4[["parm"]], CI.method = CI.methods, na.rm = na.rm, round=4)
        temp.eval[["4. imp.norm"]] <- CI.evaluator(temp.d4[["est.MI"]][["imp.norm"]], param = temp.d4[["parm"]], CI.method = CI.methods, na.rm = na.rm, round=4)
        temp.eval[["5. Bootstrap"]] <- CI.evaluator(temp.d4[["est.Dir"]], param = temp.d4[["parm"]], CI.method = Dir.methods, na.rm = na.rm, round=4)
        } #imputation==TRUE
        temp.d3[[h]] <- temp.d4 
        temp.d3[[h]]$eval <- temp.eval
        temp.d3[[h]]$sim <- ""                          ##removing raw data
        if (imputation==TRUE) {temp.d3[[h]]$MI <- ""}   ##removing MI data
      }
      saveRDS(temp.d3, paste0("sim_data","-",i,j,"-",format(Sys.time(), "%b%d"),".rds"))
      temp.d2[[j]] <- temp.d3
      rm(temp.d3)
    }
    temp.d1[[i]] <- temp.d2
    rm(temp.d2)
  }
  sim.data = temp.d1
  rm(temp.d1)
}
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
