param.finder = function(param = "beta1", start, target, digit = 1, target.digit=5, tolerance = 0.000005, seed=100, round2.nsample=40, final.nsample=400) {
  # target.digit: the number of decimal placese of the final value to be returned
  a <- Sys.time()                           # time check
  length(gregexpr("[[:digit:]]", as.character(start))[[1]])-1 
  # number of decimal places of the start
  
  # routines
  nsample = function(digit) {
    if (digit<=2) {n.sample=1}
    else if (digit==3) {n.sample=5}
    else if (digit==4) {n.sample=10}
    else {n.sample=20}
    return(n.sample)
  }
  newval = function(param0, param, err0, err, digit) {
    newval = param - sign(err)*max(round(abs(err / (abs(err0)+abs(err)) *(param - param0)), digit),10^(-digit))
  }
  
  n=10000000; n1=1000000; n2 = floor(n/n1)  # n is the sample size, n1 is size of partitions for efficiency
  err = 1         # starting value of err
  tht.mean = 1
  assign(param, start)  # for function, beta1 usually
  param.val = start     # for recording
  iter = 1
  digit = digit
  n.sample = nsample(digit)    # start from 1 sample
  
  ################### round 3 ###################
  err = 0.1                        # specific to RND3 function
  init.sign = sign(err)
  iter = 1
  result = list()
  result$parameter = param
  # result$round2 = round2           # specific to RND3 function
  result$final = rep(NA,10)
  print(Sys.time() - a)
  pb <- txtProgressBar(min=0, max = final.nsample, char = paste0("R3-",iter), style=3)
  repeat{
    print(paste("iteration:" ,iter, ", n.sample:", final.nsample))
    print(paste(param,"=",param.val))
    print(paste(final.nsample, "samples(X10mil) for the best estimate +1/-1 ticks "))
    tht <- rep(NA,final.nsample)
    if (init.sign * sign(err) < 0 | ((param.val %in% result$final) & (iter > 4))) {break}
    assign(param, param.val)   #param.val0(next to last value(before being updated))
    set.seed(seed)
    for (j in 1:final.nsample) {
      setTxtProgressBar(pb,j)
      temp <- data.frame(disease=rep(NA,n), marker=NA)
      for (i in 1:n2) {temp[((i-1)*n1+1):(i*n1),] <- datagenerator(n=n1, alpha0=alpha0, alpha1=alpha1, beta0=beta0, beta1=beta1, beta2=beta2, beta3=beta3, sig=sig, q1=q1, q2=q2, q3=q3, q4=q4, gamma=gamma, mu.V=mu.V, Sigma=Sigma, option="VDT")} ; Sys.time()-a
      tht[j] <- AUC(temp$marker[temp$disease==0], temp$marker[temp$disease==1])
      if (j %% floor(final.nsample/20) == 0) {print(mean(tht, na.rm=TRUE))}
    }
    print(tht.mean <- mean(tht)); print(range(tht)); print(tht)
    
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
param.finder(start=0.8320, target=0.8, digit=4, round2.nsample=20, final.nsample=40) #alpha0=1.6111
th.80.p70 = param.finder(start=0.83197, target=0.8, digit=4, round2.nsample=20, final.nsample=400) #alpha0=1.6111

