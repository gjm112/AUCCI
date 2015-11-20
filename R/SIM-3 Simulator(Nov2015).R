### Part II. Simulation

####################
# Structure has changed since Nov 5.
# Due to memory problem, For each single loop estimation is done. Don't wait until a large dataset is done.

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
              gamma = data.frame(q1 = c(.7, .8), q2 = c(.8, .9), q3=c(.7, .8), q4=c(.8, .9), gamma = c(.7, .8), rho = c(.5, .7)))
# alpha0 = 0 for phi==0, alpha0 = 1.6111 for phi==.7
# q3= q1, q4= q2: MAR settings

## 2.3.1.3 Param2 - Simulation ###########################################################
n = c(200, 100, 50); n.max = max(n)
n.sim = 3    # To be changed to 10,000
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
n.sim=1000
{ 
  set.seed(200)
  bgn <- Sys.time()
  # setting dimensions and seed
  d1 = 1:length(param1$alphabet$theta); d2 = 1:dim(param1$gamma)[1]; d3 = 1:length(n)
  ## steps
  imputation=TRUE    # run MI and estimates?
  Dir=FALSE          # run Direct approaches?
  #save.Est=FALSE    # save Inferences(Estimations)?
  
  
  # d1 <- 2; d2 <-1; d3 <-2; imputation=FALSE
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
  error.index =0
  
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
      rho = param1$gamma$rho[j]
        
      for (h in d3) {  #h:n=200,100,50
        pb.text <- paste0("Data+",ifelse(imputation,"MI+",""),ifelse(Dir,"Dir+",""),"Est|")
        pb <- txtProgressBar(min=0, max = n.sim, char = paste0(((i-1)*2+j-1)*3+h, "/",d123, pb.text), style=3)
        
        ## Empty shells that is outer than k (1:n.sim) 
        {
          # for A-Est: Complete data
          temp.est <- as.data.frame(matrix(NA,n.sim,(length(CI.methods)*2+1)))
          names(temp.est) <- c("AUC.hat",paste0(rep(CI.methods,each=2),c(".lb",".ub")))
          
          # for B-Est: Imputations
          temp.estMI <- list()
          for (l in MI.methods[,"methods"]) {
            temp.estMI[[l]] <- as.data.frame(matrix(NA,n.sim,(length(CI.methods)*2+1)))
            names(temp.estMI[[l]]) <- c("AUC.hat",paste0(rep(CI.methods,each=2),c(".lb",".ub")))
          }
          
          # for C-Dir: Direct Approaches
          temp.estDir.BCA <- temp.estDir.Wald <-  as.data.frame(matrix(NA,n.sim,(length(Dir.methods)*2+1)))
          names(temp.estDir.BCA) <- names(temp.estDir.Wald) <- c("AUC.hat",paste0(rep(Dir.methods,each=2),c(".lb",".ub")))
          
        }

        
        for (k in 1:n.sim) {
          # Progress bar
          setTxtProgressBar(pb,k)
          
          # A-Gen. Data generation(from 2.)
          # temp <- datagenerator.normal(n=n[h], theta=theta, phi=phi, sig.y=3)
          temp <- datagenerator(n=n[h], alpha0=alpha0, alpha1=alpha1, beta0=beta0, beta1=beta1, beta2=beta2, beta3=beta3, sig=sig, q1=q1, q2=q2, q3=q3, q4=q4, gamma=gamma, mu.V=mu.V, Sigma=Sigma, option="VDTR")
          
          # A-Est: Inferences (complete datasets)
          # pb <- txtProgressBar(min=0, max = n.sim, char = paste0(((i-1)*2+j-1)*3+h, "/",d123," (Est.com-)|"), style=3)  
          # setTxtProgressBar(pb,k)
          temp.est[k,] <- CI.i(temp,fun=AUCCI, CI.method=CI.methods, alpha=alpha, type="landscape2")
          
          # B. Imputation
          if (imputation==TRUE) {
            # B-Imp
            temp.MI <- list()
            nonimputable <- list()
            for (l in MI.methods[,"methods"]) {
              # pb <- txtProgressBar(min=0, max = n.sim, char = paste0(((i-1)*2+j-1)*3+h, "/",d123," (MI-",l,")|"), style=3)
              if (l == "pmm" |l == "logreg") { 
                temp.imp <- mice2(data=temp, method=l, m=m, predictorMatrix=cbind(0,(1 - diag(1, ncol(temp)))[,-1]), printFlag=FALSE)
                if (identical(temp.imp,"error")) {
                  temp.MI[[l]] = vector("list", m)
                  nonimputable[[l]] = TRUE
                } else {
                  temp.comp <- list()
                  for (q in 1:m) {temp.comp[[q]] = complete(temp.imp, action = q)[,c("diseaseR", "marker")]}
                  temp.MI[[l]] = temp.comp  
                  nonimputable[[l]] = FALSE
                }
              }
              else if (l == "imp.norm") {
                ###TBD!!!!!!
                temp.MI[[l]] = temp.MI[[1]]  # proxy (pmm is copied)
                nonimputable[[l]] = nonimputable[[1]]
              }
              else {print("Wrong MI method!")}
            }
            
            # B-Est (Inference on incomplete datasets with MI)
            for (l in MI.methods[,"methods"]) {
              # pb <- txtProgressBar(min=0, max = n.sim*3, char = paste0(((i-1)*2+j-1)*3+h, "/",d123," (Est.MI-",l,")|"), style=3)
              # l.index = which(l ==MI.methods[,"methods"] )
              # setTxtProgressBar(pb,k + (l.index-1)*n.sim)
              if (nonimputable[[l]] == TRUE) {
                temp.estMI[[l]][k,] <- NA
              } else {
                temp.estMI[[l]][k,] <- CI.i(temp.MI[[l]], fun=AUCCI.MI, CI.method=CI.methods, m=m, alpha=alpha, type="landscape2")
              }
              
              if (l == "pmm" & (is.na(temp.estMI[[l]][k,2])|(temp.estMI[[l]][k,2]>=1 |temp.estMI[[l]][k,2] < 0  ))) {
                error.index <- error.index +1
                temp.comp.error[[error.index]] <- temp.MI[[l]]
              }
            }
          } #imputation==TRUE
          
          # C. Direct Approach (Inferences on incomplete datasets with Direct methods(Bootstrapping))
          if (Dir==TRUE) {
            # pb <- txtProgressBar(min=0, max = n.sim, char = paste0(((i-1)*2+j-1)*3+h, "/",d123," (Est.Dir-)|"), style=3)
            # setTxtProgressBar(pb,k)
            temp.boot <- AUCCI.boot(data=temp, R=R, CI.method=Dir.methods, alpha=alpha, type="landscape2")
            temp.estDir.BCA[k,] <- cbind("-",temp.boot$CI)
            temp.estDir.Wald[k,] <- cbind("-",temp.boot$CI.Wald)
          } #Dir==TRUE
          
        }
        
        # storing A ~ C-2
        
        temp.d4 <- list(parm = data.frame(theta=theta, phi=phi, gamma=gamma, rho=rho), est.com = temp.est)
        if (imputation==TRUE) {
          temp.d4[["est.MI"]] = temp.estMI
        }
        if (Dir==TRUE)        {
          temp.d4[["est.Dir.BCA"]] = temp.estDir.BCA
          temp.d4[["est.Dir.Wald"]] = temp.estDir.Wald
        }
        
        # D-1. Evaluation of C-1(est.com)
        temp.eval <- list()
        temp.eval[["1. complete"]] <- CI.evaluator(temp.d4[["est.com"]], param = temp.d4[["parm"]], CI.method = CI.methods, na.rm = na.rm, round=4)
        names(temp.eval)[1] <- paste0("theta=",theta,", phi=",phi,", rho=",rho,", n=",n[h],", 1.complete")
        if (imputation==TRUE) {
          temp.eval[["2. pmm"]] <- CI.evaluator(temp.d4[["est.MI"]][["pmm"]], param = temp.d4[["parm"]], CI.method = CI.methods, na.rm = na.rm, round=4)
          temp.eval[["3. logreg"]] <- CI.evaluator(temp.d4[["est.MI"]][["logreg"]], param = temp.d4[["parm"]], CI.method = CI.methods, na.rm = na.rm, round=4)
          temp.eval[["4. imp.norm"]] <- CI.evaluator(temp.d4[["est.MI"]][["imp.norm"]], param = temp.d4[["parm"]], CI.method = CI.methods, na.rm = na.rm, round=4)
          names(temp.eval)[2:4] <- paste0("theta=",theta,", phi=",phi,", rho=",rho,", n=",n[h],c(", 2.pmm", ", 3.logreg", ", 4.imp.norm"))
        } #imputation==TRUE
        if (Dir==TRUE) {
          temp.eval[["5. Bootstrap-BCA"]] <- CI.evaluator(temp.d4[["est.Dir.BCA"]], param = temp.d4[["parm"]], CI.method = Dir.methods, na.rm = na.rm, round=4)
          temp.eval[["6. Bootstrap-Wald"]] <- CI.evaluator(temp.d4[["est.Dir.Wald"]], param = temp.d4[["parm"]], CI.method = Dir.methods, na.rm = na.rm, round=4)
          names(temp.eval)[ifelse(imputation==TRUE,5:6,2:3)] <- paste0("theta=",theta,", phi=",phi,", rho=",rho,", n=",n[h],c(", 5. Bootstrap-BCA", ", 6. Bootstrap-Wald"))
        } #Dir==TRUE
        temp.d3[[h]] <- temp.d4
        temp.d3[[h]]$eval <- temp.eval
      
        # Time stat
        elapsed = Sys.time()-bgn
        expected = bgn + elapsed/(((i-1)*2+j)*3+h)*(d123)
        print(paste0("bgn: ", format(bgn,"%m/%d %H:%M"), ", elapsed: ", round(elapsed,1), " min's, expected: ", format(expected,"%m/%d %H:%M"), ", i: ", paste(i,"in", length(d1)), ", j: ", paste(j,"in", length(d2)), ", h: ", paste(h,"in", length(d3)) ))
      
      } # h in d3
      saveRDS(temp.d3, paste0("sim_data","-",format(Sys.time(), "%b%d"),"-",i,j,".rds"))
      temp.d2[[j]] <- temp.d3
      
    } # j in d2
    temp.d1[[i]] <- temp.d2
  } # i in d1
  sim.data = temp.d1
  rm(temp.d1,temp.d2)
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
