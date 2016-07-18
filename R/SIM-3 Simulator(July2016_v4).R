### Part II. Simulation

########################################
library(foreach) # foreach
library(doSNOW)
library(tcltk)

exp.packages = c("foreach","MASS", "mice", "norm", "rootSolve", "sn", "tcltk")
exp.functions = names(which(eapply(globalenv(), class) == "function"))
exp.functions2 = c("Rubin", "S.stat")
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
param1 = list(alphabet = data.frame(phi = rep(c(.5,.7), each=4), 
                                    theta = c(0.7999914248, 0.89999851, 0.950002331, 0.9900000165, 0.800002545, 0.899999218, 0.950000849, 0.989999545), 
                                    theta.SE =c(5.80868E-06, 4.09076E-06, 2.76121E-06, 1.0662E-06, 7.79383E-06, 5.34159E-06, 3.42365E-06, 1.15113E-06),
                                    alpha0 = rep(c(0,1.6111), each=4), 
                                    beta1=c(0.8089, 1.4486, 1.97674, 2.96704, .8319, 1.47286, 2.00192, 2.9939)), 
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

n.sim=10000
bgn <- Sys.time()

# setting dimensions and seed
d1 = 1:length(param1$alphabet$theta); d2 = 1:dim(param1$gamma)[1]; d3 = 1:length(n)
#d1 = c(4,8)
d123 <- length(d1)*length(d2)*length(d3)

## steps
imputation=TRUE    # run MI and estimates?
Dir=FALSE          # run Direct approaches?
#save.Est=FALSE    # save Inferences(Estimations)?
na.rm = TRUE       # for evaluation
parallel = TRUE    # %dopar% or %do%?
# if parallel, should change %do% -> %dopar%, and "<<-" -> "<-" (exceptions: debug.*, count)

## 2.3.2 Simulator #######################################################################

## 2.3.2.1 sim.by.k functionalized part under i, j, h
debug.ijh <- rep(NA,3)
debug.l <- debug.k <- NA

sim.by.ijh <- function(i, j, h, d123, pb, pb.text, param, mu.V, Sigma, CI.method, alpha, n,
                       imputation, MI.method, Dir, m, Dir.method, R) {
  debug.ijh <<- c(i,j,h)
  
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
  temp.est <- as.data.frame(matrix(NA,n.sim,(length(CI.methods)*2+1)))
  names(temp.est) <- c("AUC.hat",paste0(rep(CI.methods,each=2),c(".lb",".ub")))
  
  # for B-Est: Imputations
  temp.estMI <- lapply(MI.method, function(l) {a <- as.data.frame(matrix(NA,n.sim,(length(CI.methods)*2+1))) 
                                               names(a) <- c("AUC.hat",paste0(rep(CI.methods,each=2),c(".lb",".ub"))); return(a)})
  names(temp.estMI) <- MI.method
  
  # for C-Dir: Direct Approaches
  temp.estDir.BCA <- temp.estDir.Wald <-  as.data.frame(matrix(NA,n.sim,(length(Dir.methods)*2+1)))
  names(temp.estDir.BCA) <- names(temp.estDir.Wald) <- c("AUC.hat",paste0(rep(Dir.methods,each=2),c(".lb",".ub")))
  
  # for progressbar
  
  
  #CI.method=CI.method, alpha=alpha,
  #imputation=imputation, MI.method=MI.method, m=m, Dir=Dir, Dir.method=Dir.method, R=R)
  
  # Estimation: simulation by k
  for (k in 1:n.sim) {
    debug.k <<- k
    # A-Gen. Data generation(from 2.)
    temp <- datagenerator(n=n[h], alpha0=alpha0, alpha1=alpha1, beta0=beta0, beta1=beta1, beta2=beta2, beta3=beta3, sig=sig, q1=q1, q2=q2, q3=q3, q4=q4, gamma=gamma, mu.V=mu.V, Sigma=Sigma, option="VDTR")
    # A-Est: Inferences (complete datasets)
    temp.est[k,] <- CI.i(temp,fun=AUCCI, CI.method=CI.method, alpha=alpha, type="landscape2")
    
    ####### counter ######
    count <<- count + 1
    #pb.text <- paste0("Data+",ifelse(imputation,"MI+",""),ifelse(Dir,"Dir+",""),"Est|")
    #pb <- txtProgressBar(min=0, max = n.sim, char = paste0(((i-1)*2+j-1)*3+h, "/",d123, pb.text), style=3)
    #setTxtProgressBar(pb,count)
    if (k %% 300 == 0) {write.csv(NA, paste0("progress/progress(",paste0(c(i,j,h,k),collapse="-"),")",format(Sys.time(), "%b%d-%H%M%S"),".csv"))}
    if (count %% 1000 == 0) {
      progress = count / (d123*n.sim)
      elapsed = Sys.time() - bgn
      expected = bgn + elapsed/progress
      write.csv(NA, paste0("progress/Expected(as of ",format(Sys.time(), "%b%d-%H%M"),") is ",format(expected, "%b%d-%H%M"),".csv"))
    }
    #Sys.sleep(0.01)
    #flush.console()
    ####### counter ######
    
    # B. Imputation
    if (imputation) {
      # B-Imp
      temp.MI <- list();  
      nonimputable <- lapply(MI.method,function(x) FALSE); names(nonimputable) <- MI.method
      for (l in MI.method) {
        debug.l <<- c(l,"temp.MI")
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
        debug.l <<- c(l,"temp.estMI")
        if (nonimputable[[l]]) {
          temp.estMI[[l]][k,] <- NA
        } else {
          temp.estMI[[l]][k,] <- CI.i(temp.MI[[l]], fun=AUCCI.MI, CI.method=CI.method, m=m, alpha=alpha, type="landscape2")
        }
      }
    } #imputation==TRUE
    debug.l <<- c(l,NA)
    
    # C. Direct Approach (Inferences on incomplete datasets with Direct methods(Bootstrapping))
    if (Dir) {
      temp.boot <- AUCCI.boot(data=temp, R=R, CI.method=Dir.methods, alpha=alpha, type="landscape2")
      temp.estDir.BCA[k,] <- cbind("-",temp.boot$CI)
      temp.estDir.Wald[k,] <- cbind("-",temp.boot$CI.Wald)
    } #Dir==TRUE
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
count <- 1   # for time checking in dopar

pb <- txtProgressBar(min=0, max = length(d1)*length(d2), style=3)
for (i in d1) {
#  if (i > 6) {    # break and resume by controlling numbers
  for (j in d2) {
    if (parallel) {cl <- makeSOCKcluster(4)}      # for parallel
    if (parallel) {registerDoSNOW(cl)}            # for parallel
    print(c(i,j, paste0("out of ", length(d1),"x", length(d2))))
    sim.data  <- temp.d1 <- temp.d2 <- list()
    bgn <- Sys.time()
    # progress(h)
    setTxtProgressBar(pb,(i-1)*length(d2)+j)
    
    if (parallel) {
      temp.d3  <- foreach (h=d3, .packages=exp.packages, .export=exp.functions2) %dopar% {
        set.seed(i*j*h*100)
        sim.by.ijh (i=i, j=j, h=h, d123=d123, pb=pb, pb.text=pb.text, param=param1, mu.V=mu.V, Sigma=Sigma, CI.method=CI.methods, alpha=alpha, n=n,
                    imputation=imputation, MI.method = MI.methods[,"methods"], 
                    Dir=Dir, m=m, Dir.method=Dir.methods, R=R)
      }      
    } else {
      temp.d3  <- foreach (h=d3, .packages=exp.packages, .export=exp.functions2) %do% {
        set.seed(i*j*h*100)
        sim.by.ijh (i=i, j=j, h=h, d123=d123, pb=pb, pb.text=pb.text, param=param1, mu.V=mu.V, Sigma=Sigma, CI.method=CI.methods, alpha=alpha, n=n,
                    imputation=imputation, MI.method = MI.methods[,"methods"], 
                    Dir=Dir, m=m, Dir.method=Dir.methods, R=R)
      }
    }

    if (parallel) {stopCluster(cl)}               # for parallel
    saveRDS(temp.d3, paste0("R/Simdata3/sim_data","-",format(Sys.time(), "%b%d"),"-",i,j,".rds"))
    # Time stat
    elapsed = Sys.time()-bgn
    expected = bgn + elapsed/((i-1)*2+j)*12
    print(paste0("bgn: ", format(bgn,"%m/%d %H:%M"), ", elapsed: ", round(elapsed,1), " min's, expected: ", format(expected,"%m/%d %H:%M"), ", i: ", paste(i,"in", length(d1)), ", j: ", paste(j,"in", length(d2)) ))
    
    temp.d2[[j]] <- temp.d3
  } # j in d2
  temp.d1[[i]] <- temp.d2
# } # break and resume
} # i in d1
# sim.data = temp.d1; rm(temp.d1,temp.d2)

