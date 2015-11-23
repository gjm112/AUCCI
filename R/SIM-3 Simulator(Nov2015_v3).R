### Part II. Simulation

########################################
library(foreach) # foreach
library(doMC)    # %dopar% in foreach
registerDoMC(cores=4)
library(doSNOW)
library(tcltk)
cl <- makeSOCKcluster(3)
registerDoSNOW(cl)

exp.packages = c("foreach","MASS", "mice", "norm", "rootSolve", "sn", "tcltk")
exp.functions = names(which(eapply(globalenv(), class) == "function"))
########################################

# Structure has changed since Nov 5.
# Due to memory problem, For each single loop estimation is done. Don't wait until a large dataset is done.

## 2.3.1. Parameters  ###################################################################
## 2.3.1.1 Param1 - Distributional #######################################################
# parameters (except alpha0 and beta1)
mu.V = rep(0,5)
Sigma = matrix(c(1, 0, 0.3, 0.4, -0.4, 0, 1, 0.2, 0.2, 0, 0.3, 0.2, 1, 0.7, -0.5, 0.4, 0.2, .7, 1, -0.2, -0.4, 0, -0.5, -0.2, 1), 5,5)
alpha1 = rep(1,5)
beta0 = 0
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

## 2.3.1.2 Param2 - Simulation ###########################################################
n = c(200, 100, 50)
alpha = 0.05
# CI.methods from 1-1 [CI] Base functions.R
CI.methods = c("Bm", "HM1", "HM2", "NS1", "NW", "NS2", "Mee", "DL", "RG", "DB", "DG", "CM")
# MI.methods from 1-5 [CI] AUCCI_MI
MI.methods = data.frame(functions=c(rep("mice2",2),rep("MI.norm",3)), methods=c("pmm","logreg","simple","coinflip","adaptive"))
m = 10
# Direct methods with bootstrap. R = # of resamples
Dir.methods = c("naive", "BG", "MS", "SP", "IPW", "He")
R = 30

n.sim=10
bgn <- Sys.time()
# setting dimensions and seed
d1 = 1:length(param1$alphabet$theta); d2 = 1:dim(param1$gamma)[1]; d3 = 1:length(n)
d123 <- length(d1)*length(d2)*length(d3)

## steps
imputation=TRUE    # run MI and estimates?
Dir=FALSE          # run Direct approaches?
#save.Est=FALSE    # save Inferences(Estimations)?
na.rm = TRUE       # for evaluation

## 2.3.2 Simulator #######################################################################

## 2.3.2.1 sim.by.k functionalized part under i, j, h

# Those below are global variables, and input and output variables: 
# They are the final results but are not explicitly returned.
# High level assigning operator (<<-) is used.
# temp.est, temp.estMI, temp.estDir.BCA, temp.estDir.Wald
sim.by.k = function(k, h, n,
                    alpha0, alpha1, beta0, beta1, beta2, beta3, sig, q1, q2, q3, q4, gamma, rho, mu.V, Sigma, 
                    CI.method, alpha,
                    imputation, MI.method, m, Dir, Dir.method, R) {
  # A-Gen. Data generation(from 2.)
  temp <- datagenerator(n=n[h], alpha0=alpha0, alpha1=alpha1, beta0=beta0, beta1=beta1, beta2=beta2, beta3=beta3, sig=sig, q1=q1, q2=q2, q3=q3, q4=q4, gamma=gamma, mu.V=mu.V, Sigma=Sigma, option="VDTR")
  # A-Est: Inferences (complete datasets)
  temp.est[k,] <<- CI.i(temp,fun=AUCCI, CI.method=CI.method, alpha=alpha, type="landscape2")
  
  # B. Imputation
  if (imputation) {
    # B-Imp
    temp.MI <- list();  
    nonimputable <- lapply(MI.method,function(x) FALSE); names(nonimputable) <- MI.method
    for (l in MI.method) {
      if (l == "pmm" |l == "logreg") { 
        temp.MI[[l]] <- mice2(data=temp[,-1], method=l, m=m, printFlag=FALSE)
        if (is.na(temp.MI[[l]][[1]][[1]][[1]])) { nonimputable[[l]] = TRUE}
      }
      else if (l == "simple"|l == "coinflip"|l == "adaptive") {
        temp.MI[[l]] <- MI.norm(data=temp[,-1], rounding=l, m=m, showits=FALSE)
      }
      else {stop("Wrong MI method!")}
    }
    
    # B-Est (Inference on incomplete datasets with MI)
    for (l in MI.method) {
      if (nonimputable[[l]]) {
        temp.estMI[[l]][k,] <<- NA
      } else {
        temp.estMI[[l]][k,] <<- CI.i(temp.MI[[l]], fun=AUCCI.MI, CI.method=CI.method, m=m, alpha=alpha, type="landscape2")
      }
    }
  } #imputation==TRUE
  
  # C. Direct Approach (Inferences on incomplete datasets with Direct methods(Bootstrapping))
  if (Dir) {
    temp.boot <- AUCCI.boot(data=temp, R=R, CI.method=Dir.methods, alpha=alpha, type="landscape2")
    temp.estDir.BCA[k,] <<- cbind("-",temp.boot$CI)
    temp.estDir.Wald[k,] <<- cbind("-",temp.boot$CI.Wald)
  } #Dir==TRUE
}

sim.by.ijh <- function(i, j, h, d123, param, mu.V, Sigma, CI.method, alpha, n,
                       imputation, MI.method, Dir, m, Dir.method, R) {
  
  ### 1. basic setting
  # param(i)
  alpha0 = param$alphabet$alpha0[i]
  beta1  = param$alphabet$beta1[i]
  phi    = param$alphabet$phi[i]
  theta  = param$alphabet$theta[i]
  # param(j)
  gamma  = param$gamma$gamma[j]
  q1  = param$gamma$q1[j]
  q2  = param$gamma$q2[j]
  q3  = param$gamma$q3[j]
  q4  = param$gamma$q4[j]
  rho = param$gamma$rho[j]
  
  ### 2. Empty shells that is outer than k (1:n.sim) 
  # for A-Est: Complete data
  temp.est <<- as.data.frame(matrix(NA,n.sim,(length(CI.methods)*2+1)))
  names(temp.est) <<- c("AUC.hat",paste0(rep(CI.methods,each=2),c(".lb",".ub")))
  
  # for B-Est: Imputations
  temp.estMI <<- lapply(MI.method, function(l) {a <- as.data.frame(matrix(NA,n.sim,(length(CI.methods)*2+1))) 
                                               names(a) <- c("AUC.hat",paste0(rep(CI.methods,each=2),c(".lb",".ub"))); return(a)})
  names(temp.estMI) <<- MI.method
  
  # for C-Dir: Direct Approaches
  temp.estDir.BCA <<- temp.estDir.Wald <<-  as.data.frame(matrix(NA,n.sim,(length(Dir.methods)*2+1)))
  names(temp.estDir.BCA) <<- names(temp.estDir.Wald) <<- c("AUC.hat",paste0(rep(Dir.methods,each=2),c(".lb",".ub")))
  
  # Estimation: simulation by k
  b <- foreach (k=1:n.sim) %do% {  # can't do dopar. dopar doesn't produce values for temp.est etc
    # setTxtProgressBar(pb,k, title=pb.text)
    # sim.by.k function is manipulating the data by using High level assigning operator (<<-): 
    # temp.est, temp.estMI, temp.estDir.BCA, temp.estDir.Wald
    sim.by.k (k=k, h=h, n=n,
              alpha0=alpha0, alpha1=alpha1, beta0=beta0, beta1=beta1, beta2=beta2, beta3=beta3, 
              sig=sig, q1=q1, q2=q2, q3=q3, q4=q4, gamma=gamma, rho=rho, mu.V=mu.V, Sigma=Sigma,
              CI.method=CI.methods, alpha=alpha,
              imputation=imputation, MI.method=MI.method, m=m, Dir=Dir, Dir.method=Dir.methods, R=R)
  }
  
  # storing A ~ C-2
  temp.d4 <- list(parm = data.frame(theta=theta, phi=phi, gamma=gamma, rho=rho, n=n[h]), est.com = temp.est)
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
  names(temp.eval)[1] <- paste0("theta=",round(theta,2),", phi=",phi,", rho=",rho,", n=",n[h],", 1.complete")
  if (imputation==TRUE) {
    temp.eval <- c(temp.eval, lapply(temp.d4[["est.MI"]], function(x) CI.evaluator(x, param = temp.d4[["parm"]], CI.method = CI.methods, na.rm = na.rm, round=4) ))
    names(temp.eval)[2:6] <- paste0("theta=",round(theta,2),", phi=",phi,", rho=",rho,", n=",n[h],c(", 2.pmm", ", 3.logreg", ", 4.simple", ", 5.coinflip", ", 6.adaptive"))
  } #imputation==TRUE
  if (Dir==TRUE) {
    temp.eval[["7. Bootstrap-BCA"]] <- CI.evaluator(temp.d4[["est.Dir.BCA"]], param = temp.d4[["parm"]], CI.method = Dir.methods, na.rm = na.rm, round=4)
    temp.eval[["8. Bootstrap-Wald"]] <- CI.evaluator(temp.d4[["est.Dir.Wald"]], param = temp.d4[["parm"]], CI.method = Dir.methods, na.rm = na.rm, round=4)
    names(temp.eval)[if(imputation) {c(7,8)} else{c(2,3)}] <- paste0("theta=",round(theta,2),", phi=",phi,", rho=",rho,", n=",n[h],c(", 7. Bootstrap-BCA", ", 8. Bootstrap-Wald"))
  } #Dir==TRUE
  temp.d4$eval <- temp.eval
  return(temp.d4)
}


## 2.3.2.2 Simulation package (Data generation + Inference + Evaluation)  ################

{ 
  # data structure
  # temp.d1 (=sim.data) = list of 6 elements(temp.d2): one for each theta & phi
  # temp.d2             = list of 2 elements(unnamed list): one for each rho
  # temp.d3             = list of 4 elements(temp.sim, param, temp.est, temp.eval)
  # temp.sim            = list of 10000(=n.sim) elements: the unit samples
  # temp.est            = list of 3 elements(temp.est.n): one for n(200, 100, 5)
  # temp.est.n          = a dataframe of estimators(point and CI's)
  # temp.eval           = a dataframe of evaluation of the estimators
  
  # creating empty lists : remove when final for loop is developed and they are not necessary any more.
  sim.data  <- temp.d1 <- temp.d2 <- list()
  i <- j <- 1
  a <- Sys.time()
  temp.d3  <- foreach (h=d3, .packages=exp.packages, .export=exp.functions) %dopar% {
    set.seed(i*j*h*100)
    # progress(h)
    # pb.text <- paste0("Data+",ifelse(imputation,"MI+",""),ifelse(Dir,"Dir+",""),"Est|")
    # pb <- txtProgressBar(min=0, max = n.sim, char = paste0(((i-1)*2+j-1)*3+h, "/",d123, pb.text), style=3)
    # sim.by.k function is manipulating the data by using High level assigning operator (<<-): 
    # temp.est, temp.estMI, temp.estDir.BCA, temp.estDir.Wald
    sim.by.ijh (i=i, j=j, h=h, d123=d123, param=param1, mu.V=mu.V, Sigma=Sigma, CI.method=CI.methods, alpha=alpha, n=n,
                imputation=imputation, MI.method = MI.methods[,"methods"], 
                Dir=Dir, m=m, Dir.method=Dir.methods, R=R)
  }
  Sys.time() - a
  saveRDS(temp.d3, paste0("sim_data","-",format(Sys.time(), "%b%d"),"-",i,j,".rds"))
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

