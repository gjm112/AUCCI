### Part II. Simulation
## 2-2 [SIM] Evaluator: 2.2.1 CI.evaluator ################################################

## 2.2.1 CI.evaluator #####################################################################
CI.evaluator = function(data, param, CI.method = sub(".lb","",matrix(names(data)[-1],nrow=2)[1,]),
                        na.rm = TRUE, round=Inf) {
  # data: a data.frame of CI.method by n.sim
  # param: a data.frame that has theta values
  # round: # of decimal points
  
  msmt = c("CP","LNCP","RNCP","CIL", "ZWI","NaN")
  
  # data formatting
  len.col = length(CI.method); len.row = length(msmt)
  eval = as.data.frame(matrix(NA,len.row,len.col))
  names(eval) = CI.method
  rownames(eval) = msmt
  
  # informations
  theta = param$theta
  
  # CP calculation
  temp.1 <- temp <- data[,-1];
  temp.1[,] <- 0
  # 1 if lowerbound captures theta, 0 otherwise
  temp.1[,seq(1,2*len.col-1,by=2)] = (temp[,seq(1,2*len.col-1,by=2)] < theta)
  # 1 if upperbound captures theta, 0 otherwise
  temp.1[,seq(2,2*len.col,by=2)] = (temp[,seq(2,2*len.col,by=2)] >= theta)
  
  temp.2 = temp.1[,seq(1,2*len.col-1,by=2)]
  names(temp.2) = CI.method
  
  for (i in 1:len.col) {temp.2[,i] = temp.1[,2*i-1]*temp.1[,2*i]}
  # CP
  eval[1,] = apply(temp.2,2,mean,na.rm=na.rm)
  # NCP
  eval[2:3,] = matrix(1 - apply(temp.1, 2, mean, na.rm=na.rm),nrow=2)
  # CIL
  temp.3 = temp[,seq(2,2*len.col,by=2)] - temp[,seq(1,2*len.col-1,by=2)]
  eval[4,] = apply(temp.3, 2, mean, na.rm=na.rm)
  # ZWI
  eval[5,] = apply((temp.3 == 0), 2, mean, na.rm=na.rm)
  # NaN prob.
  temp.4 = is.na(temp.1[,seq(1,len.col*2-1,by=2)])*is.na(temp.1[,seq(2,len.col*2,by=2)])
  eval[6,] = apply(temp.4, 2, mean)
  
  return(round(eval,round))
  # if want to remove scientific notations,
  # apply "format(xx, scientific=F)" to each single cells, not to the whole dataframe(subject to error)
  # or manipulate the global setting: option(scipen = 999) and then return it by option(scipen = 0).
}