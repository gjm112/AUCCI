param.finder = function(param = "beta1", start, target, digit = 1, target.digit=5, tolerance = 0.000005, seed=100) {
  # target.digit: the number of decimal placese of the final value to be returned
  a <- Sys.time()                           # time check
  length(gregexpr("[[:digit:]]", as.character(start))[[1]])-1 
                                            # number of decimal places of the start
  
  # routines
  nsample = function(digit) {
    if (digit<=2) {n.sample=1}
    else if (digit==3) {n.sample=5}
    else if (digit==4) {n.sample=10}
    else {n.sample=50}
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
  itin = 1
  digit = digit
  n.sample = nsample(digit)    # start from 1 sample
  
  
  repeat{
    print(paste("iteration:" ,itin, ", n.sample:", n.sample))
    if (err < tolerance & digit == target.digit+1) {break}
    
    tht = rep(NA,n.sample)
    set.seed(seed)
    print(paste(param,"=",param.val))
    for (j in 1:n.sample) {
      temp <- data.frame(disease=rep(NA,n), marker=NA)
      for (i in 1:n2) {temp[((i-1)*n1+1):(i*n1),] <- datagenerator(n=n1, alpha0=alpha0, alpha1=alpha1, beta0=beta0, beta1=beta1, beta2=beta2, beta3=beta3, sig=sig, q1=q1, q2=q2, q3=q3, q4=q4, gamma=gamma, mu.V=mu.V, Sigma=Sigma, DT.only=TRUE)} ; Sys.time()-a
      tht[j] <- AUC(temp$marker[temp$disease==0], temp$marker[temp$disease==1]); Sys.time()-a
    }
    print(tht.mean <- mean(tht)); print(range(tht)); print(tht)
    err0 = err                # old value
    err = tht.mean - target   # new value
    param.val0 = param.val    # old value
    if (itin != 1 & err0 * err < 0 ) {    # if changes sign, then increase digit
      digit = digit + 1
      n.sample = nsample(digit)
      param.val = newval(param0=param.val0, param=param.val, err0=err0, err=err, digit=digit)
      assign(param, param.val)    # new value
    } else {
      multiplier = ifelse(abs(err0/(err0-err))>=3,max(floor(abs(err0/(err0-err))/2),1), 1)
      param.val = param.val - multiplier*sign(err)*10^(-digit)
      assign(param, param.val)    # new value
    }
    itin = itin + 1
  }
  
  print(Sys.time() - a)
  
  ################### round 2 ###################
  
  print("500 samples(X10mil) for +/- param estimate")
  ### 500 samples (param +/- 1digit)
  ## 500 samples (param +0)
  assign(param, param.val)
  set.seed(seed)
  for (j in 1:500) {
    temp <- data.frame(disease=rep(NA,n), marker=NA)
    for (i in 1:n2) {temp[((i-1)*n1+1):(i*n1),] <- datagenerator(n=n1, alpha0=alpha0, alpha1=alpha1, beta0=beta0, beta1=beta1, beta2=beta2, beta3=beta3, sig=sig, q1=q1, q2=q2, q3=q3, q4=q4, gamma=gamma, mu.V=mu.V, Sigma=Sigma, DT.only=TRUE)} ; Sys.time()-a
    tht[j] <- AUC(temp$marker[temp$disease==0], temp$marker[temp$disease==1])
  }
  print(paste(param,"=",param.val)); print(tht.mean <- mean(tht)); print(range(tht)); print(tht)
  
  ## 500 samples (param -1)
  assign(param, param.val-10^-5)
  set.seed(seed)
  for (j in 1:500) {
    temp <- data.frame(disease=rep(NA,n), marker=NA)
    for (i in 1:n2) {temp[((i-1)*n1+1):(i*n1),] <- datagenerator(n=n1, alpha0=alpha0, alpha1=alpha1, beta0=beta0, beta1=beta1, beta2=beta2, beta3=beta3, sig=sig, q1=q1, q2=q2, q3=q3, q4=q4, gamma=gamma, mu.V=mu.V, Sigma=Sigma, DT.only=TRUE)} ; Sys.time()-a
    tht[j] <- AUC(temp$marker[temp$disease==0], temp$marker[temp$disease==1]); Sys.time()-a
  }
  print(paste(param,"=",param.val-10^-5)); print(tht.mean <- mean(tht)); print(range(tht)); print(tht)
  
  
  ## 500 samples (param +1)
  assign(param, param.val+10^-5)
  set.seed(seed)
  for (j in 1:500) {
    temp <- data.frame(disease=rep(NA,n), marker=NA)
    for (i in 1:n2) {temp[((i-1)*n1+1):(i*n1),] <- datagenerator(n=n1, alpha0=alpha0, alpha1=alpha1, beta0=beta0, beta1=beta1, beta2=beta2, beta3=beta3, sig=sig, q1=q1, q2=q2, q3=q3, q4=q4, gamma=gamma, mu.V=mu.V, Sigma=Sigma, DT.only=TRUE)} ; Sys.time()-a
    tht[j] <- AUC(temp$marker[temp$disease==0], temp$marker[temp$disease==1]); Sys.time()-a
  }
  print(paste(param,"=",param.val+10^-5)); print(tht.mean <- mean(tht)); print(range(tht)); print(tht)
  
  print(Sys.time() - a)
}