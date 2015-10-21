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
                                    theta = c(0.7999914248, 0.89999851, 0.950002331, 0.800002545, 0.899999218, .95), 
                                    theta.SE =c(5.80868E-06,4.09076E-06,2.76121E-06,7.79383E-06, 5.34159E-06, 0) ,
                                    alpha0 = rep(c(0,1.6111), each=3), 
                                    beta1=c(0.8089, 1.4486, 1.97674,.8319,1.47286,2.0024)), 
              gamma = data.frame(q1 = c(.8, .95), q2 = c(.9, .95), q3=q1, q4=q2, gamma = c(.8, .95)))
# alpha0 = 0 for phi==0, alpha0 = 1.6111 for phi==.7
# q3= q1, q4= q2: MAR settings

## 2.3.1.3 Param2 - Simulation ###########################################################
n = c(200, 100, 50); n.max = max(n)
n.sim = 1000    # To be changed to 10,000
alpha = c(0.1, 0.05, 0.01)

# CI.methods from 1-1 [CI] Base functions.R
CI.methods = c("Bm", "HM1", "HM2", "NS1", "NW", "NS2", "Mee", "DL", "RG", "DB", "DG", "CM")

AUC.fun = AUCCI

## 2.3.2 Simulator #######################################################################

## 2.3.2.1 Simulation package (Data generation + Inference + Evaluation)  ################
{ 
  # setting dimensions and seed
  d1 = length(param1$alphabet$theta); d2 = dim(param1$gamma)[1]; d3 = length(n); set.seed = 50
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
  sim.data  <- temp.d1 <- temp.d2 <- temp.d3 <- temp.sim <- temp <- list() 
  temp.est <- list()
  temp.est.n <- as.data.frame(matrix(NA,n.sim,(length(CI.methods)*2+1)))
  names(temp.est.n) <- c("AUC.hat",paste0(rep(CI.methods,each=2),c(".lb",".ub")))
  
  # loop
  for (i in 1:d1) {
    alpha0 = param1$alphabet$alpha0[i]
    beta1  = param1$alphabet$beta1[i]
    phi    = param1$alphabet$phi[i]
    theta  = param1$alphabet$theta[i]
    
    for (j in 1:d2){
      gamma  = param1$gamma$gamma[j]
      q1  = param1$gamma$q1[j]
      q2  = param1$gamma$q2[j]
      q3  = param1$gamma$q3[j]
      q4  = param1$gamma$q4[j]
      pb <- txtProgressBar(min=0, max = n.sim, char = paste0((i-1)*2+j, "/",d1*d2," (DataGen)|"), style=3)
      for (k in 1:n.sim) {
        # A. Data generation(from 2.)
        temp.sim[[k]] <- temp <- datagenerator(n=n.max, alpha0=alpha0, alpha1=alpha1, beta0=beta0, beta1=beta1, beta2=beta2, beta3=beta3, sig=sig, q1=q1, q2=q2, q3=q3, q4=q4, gamma=gamma, mu.V=mu.V, Sigma=Sigma, option="VDTR")
        setTxtProgressBar(pb,k)
      }
      pb <- txtProgressBar(min=0, max = n.sim*length(n), char = paste0((i-1)*2+j, "/",d1*d2," (Estimate)|"), style=3)
      for (p in 1:d3) {      # 200, 100, 50
        len = n[p]
        for (k in 1:n.sim) {
          # B. Inferences(from 1.)
          temp.est.n[k,]<- CI.i(temp.sim[[k]][1:len,],fun=AUC.fun, CI.method=CI.methods, type="landscape2")
          setTxtProgressBar(pb,k + (p-1)*n.sim)
        }
        temp.est[[p]] = temp.est.n
      }

      # storing A and B
      temp.d2[[j]] <- temp.d3 <- list(sim = temp.sim, 
                                      parm = data.frame(theta=theta, phi=phi, gamma=gamma),
                                      est = temp.est) 
      # C. Evaluation of B
      temp.eval <- list()
      for (p in 1:length(n)) {      # 200, 100, 50
        temp.eval[[p]] <- CI.evaluator(temp.d3, param = 2, est = 3, n.i = p, CI.method = CI.methods, na.rm = na.rm)      
      }
      temp.d2[[j]]$eval <- temp.eval
    }
    temp.d1[[i]] <- temp.d2
  }
  sim.data = temp.d1
  rm(temp.d1)
}

## Optional: saving the datafile  ########################################################
saveRDS(sim.data, "sim_data_0927.rds")
sim.data.0927 <- readRDS("sim_data_0927.rds")

## 2.3.2.2 Simulation by part (Evaluation only)  #########################################
{
  # d1 = 6; d2 = 2;na.rm=ra.rm; set.seed=...
  # loop
  set.seed(100)
  for (i in 1:d1) {
    for (j in 1:d2){ 
      for (p in 1:d3) {
        sim.data[[i]][[j]]$eval[[p]] <- CI.evaluator(sim.data[[i]][[j]], param = 2, est = 3, n.i = p, CI.method = CI.methods, na.rm = na.rm)
      }
    }
  }
}


## Example: checking prevalence rate #####################################################
for(i in 1:20) {print(mean(sim.data[[2]][[2]][[1]][[i]]$R))}

## addressing example
# phi = .7, theta=.8  (=> 4th)
# rho = .5            (=> 1th)
# eval                (=> 4th)    # or 5th sample data (=> [[1]][[5]])
sim.data[[4]][[1]][[4]]
head(sim.data[[4]][[1]][[1]][[5]])
head(sim.data[[4]][[1]][[3]])

t(sim.data[[4]][[2]][[4]][[1]])[,1]
