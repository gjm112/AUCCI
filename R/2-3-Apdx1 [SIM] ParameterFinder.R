### Part II. Simulation
## 2-3 Apdx1 [SIM] Simulator - Appendix. parameter finder


#requires MASS, AUC, datagenerator
#########################################

param.finder = function(param = "beta1", start, target, digit = 1, target.digit=5, tolerance = 0.000005, seed=100, round2.nsample=40, final.nsample=Inf, round=1) {
  # target.digit: the number of decimal placese of the final value to be returned
  # round: starting round among 1~3
  a <- Sys.time()                           # time check
  length(gregexpr("[[:digit:]]", as.character(start))[[1]])-1 
                                            # number of decimal places of the start
  
  # routines
  nsample = function(digit, stdev=.0001) {
    if (digit<=2) {n.sample=1}
    else if (digit==3) {n.sample=min(round((10^4 * stdev)^2,0),5)}
    else if (digit==4) {n.sample=min(round((5*10^4 * stdev)^2,-1),10)}
    else { n.sample=max(min(round((stdev/tolerance)^2,-2),20), 50) }
    return(n.sample)
  }
  newval = function(param0, param, err0, err, digit) {
    newval = param - sign(err)*max(round(abs(err / (abs(err0)+abs(err)) *(param - param0)), digit),10^(-digit))
  }

  n=10000000; n1=1000000; n2 = floor(n/n1)  # n is the sample size, n1 is size of partitions for efficiency
  if (round==1) {err = 1}   # starting value of err
  else if (round==2) {err = .5}
  else {err=.1}
  init.sign = sign(err)
  tht.mean = 1
  assign(param, start)  # for function, beta1 usually
  param.val = start     # for recording
  iter = 1
  digit = digit
  n.sample = nsample(digit)    # start from 1 sample
  
  result = list()
  result$parameter = param
  
  a <- Sys.time()
  
  if (round <= 1) {
    repeat{
      print(paste("iteration:" ,iter, ", n.sample:", n.sample))
      if (err < tolerance/3 & digit == target.digit+1) {break}
      
      tht = rep(NA,n.sample)
      set.seed(seed)
      print(paste(param,"=",param.val))
      pb <- txtProgressBar(min=0, max = n.sample, char = paste0("R1"), style=3)
      for (j in 1:n.sample) {
        setTxtProgressBar(pb,j)
        temp <- data.frame(disease=rep(NA,n), marker=NA)
        for (i in 1:n2) {temp[((i-1)*n1+1):(i*n1),] <- datagenerator(n=n1, alpha0=alpha0, alpha1=alpha1, beta0=beta0, beta1=beta1, beta2=beta2, beta3=beta3, sig=sig, q1=q1, q2=q2, q3=q3, q4=q4, gamma=gamma, mu.V=mu.V, Sigma=Sigma, option="VDT")} ; Sys.time()-a
        tht[j] <- AUC(temp$marker[temp$disease==0], temp$marker[temp$disease==1]); Sys.time()-a
      }
      print(tht.mean <- mean(tht)); print(range(tht)); print(tht)
      err0 = err                # old value
      err = tht.mean - target   # new value
      stdev = sd(tht)
      param.val0 = param.val    # old value
      if (iter != 1 & err0 * err < 0 ) {    # if changes sign, then increase digit
        digit = digit + 1
        n.sample = nsample(digit, stdev)
        param.val = newval(param0=param.val0, param=param.val, err0=err0, err=err, digit=digit)
        assign(param, param.val)    # new value
      } else {
        multiplier = ifelse(abs(err0) < abs(err), 1, ifelse(abs(err0/(err0-err))>=3,min(max(floor(abs(err0/(err0-err))),1),9), 1))        
        param.val = param.val - multiplier*sign(err)*10^(-digit)
        assign(param, param.val)    # new value
      }
      iter = iter + 1
      err = tht.mean - target
    }
    
    print(Sys.time() - a)
  }

  
  ################### round 2 ###################
  if (round <= 2) {
    if (round == 2) { param.val0 = start; stdev = 0.0001}
    param.val = param.val0
    round2 <- data.frame(param = param.val+(1:7-4)*10^(-target.digit), theta.hat = NA)
    n.sample = min(max(round((stdev/tolerance/10)^2,-2),20), round2.nsample)
    print(paste(n.sample, "samples(X10mil) for the best estimate +3/-3 ticks "))
    
    for (k in 1:7){
      assign(param, round2$param[k])
      set.seed(seed)
      print(paste(param,"=",round2$param[k]))
      pb <- txtProgressBar(min=0, max = n.sample, char = paste0("R2-",k), style=3)
      for (j in 1:n.sample) {
        setTxtProgressBar(pb,j)
        temp <- data.frame(disease=rep(NA,n), marker=NA)
        for (i in 1:n2) {temp[((i-1)*n1+1):(i*n1),] <- datagenerator(n=n1, alpha0=alpha0, alpha1=alpha1, beta0=beta0, beta1=beta1, beta2=beta2, beta3=beta3, sig=sig, q1=q1, q2=q2, q3=q3, q4=q4, gamma=gamma, mu.V=mu.V, Sigma=Sigma, option="VDT")} ; Sys.time()-a
        tht[j] <- AUC(temp$marker[temp$disease==0], temp$marker[temp$disease==1])
        if (j %% floor(n.sample/5) == 0) {print(mean(tht, na.rm=TRUE))}
      }
      print(tht.mean <- mean(tht)); print(range(tht)); print(tht)
      stdev = sd(tht)
      round2$theta.hat[k] = tht.mean 
    }
    print(round2)
    lmout = lm(theta.hat ~ param, data=round2)
    param.val <- round((target - lmout$coef[1]) / lmout$coef[2], target.digit)
    print(paste("The best guess is", param.val))
    print(Sys.time() - a)   
    result$round2 = round2
  }

 
  ################### round 3 ###################
  
  if (round >= 3) { param.val0 = start; stdev = 0.0001}
  else {
    if (param.val %in% round2$param) {
      err = round2$theta.hat[which(round2$param == param.val)] - target   
    } else {err = 0.1}
  }
  init.sign = sign(err)
  iter = 1
  n.sample = min(max(round((stdev/tolerance)^2,-2),300), final.nsample)
  print(paste(n.sample, "samples(X10mil) for the best estimate +1/-1 ticks "))
  result$final = rep(NA,10)
  
  repeat{
    pb <- txtProgressBar(min=0, max = n.sample, char = paste0("R3-",iter), style=3)
    n.sample = min(max(round((stdev/tolerance)^2,-2),300), final.nsample)
    print(paste("iteration:" ,iter, ", n.sample:", n.sample))
    if (init.sign * sign(err) < 0 | ((param.val %in% result$final) & (iter > 4))) {break}
    tht <- rep(NA,n.sample)
    print(paste(param,"=",param.val))
    assign(param, param.val)   #param.val0(next to last value(before being updated))
    set.seed(seed)
    for (j in 1:n.sample) {
      setTxtProgressBar(pb,j)
      temp <- data.frame(disease=rep(NA,n), marker=NA)
      for (i in 1:n2) {temp[((i-1)*n1+1):(i*n1),] <- datagenerator(n=n1, alpha0=alpha0, alpha1=alpha1, beta0=beta0, beta1=beta1, beta2=beta2, beta3=beta3, sig=sig, q1=q1, q2=q2, q3=q3, q4=q4, gamma=gamma, mu.V=mu.V, Sigma=Sigma, option="VDT")} ; Sys.time()-a
      tht[j] <- AUC(temp$marker[temp$disease==0], temp$marker[temp$disease==1])
      if (j %% floor(n.sample/20) == 0) {print(mean(tht, na.rm=TRUE))}
    }
    print(tht.mean <- mean(tht)); print(range(tht)); print(tht)
    stdev = sd(tht)

    # storing the last values to result
    result$final[iter] = param.val  #recording beta1 values
    result[[iter+3]] = list(param.val = param.val, thetas = tht)
    # for the next iteration
    iter = iter + 1
    err = tht.mean - target
    param.val = param.val -sign(err)*10^(-target.digit)    
  }
  
  return(result)

}



# th.80 = param.finder(start=0.80889, target=0.8, digit=5)
# th.90 = param.finder(start=1.6, target=0.9, digit=1)
# th.80b = param.finder(start=0.8088, target=0.8, digit=4, round2.nsample=20, final.nsample=40)
param.finder(start=1.473, target=0.9, digit=3,round=1) #alpha0=1.6111

