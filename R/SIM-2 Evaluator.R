### Part II. Simulation
## 2-2 [SIM] Evaluator: 2.2.1 CI.evaluator ################################################

## 2.2.1 CI.evaluator #####################################################################
CI.evaluator = function(data, param, CI.method = sub(".lb","",matrix(names(data)[starting.col:length(names(data))],nrow=2)[1,]),
                        na.rm = TRUE, round=Inf, starting.col = 2) {
  # data: a data.frame of CI.method by n.sim
  # param: a data.frame that has theta values
  # round: # of decimal points
  # starting.col: to exclude 1st col (where usually point estimate is)
  
  msmt = c("CP","LNCP","RNCP","NaN","CIL", "ZWI")
  
  # data formatting
  len.col = length(CI.method); len.row = length(msmt)
  eval = as.data.frame(matrix(NA,len.row,len.col))
  names(eval) = CI.method
  rownames(eval) = msmt
  
  # informations
  theta = param$theta
  n.sim = dim(data)[1]
  
  # excluding point estimate and start from interval estimates.
  if (starting.col==1) {
    temp.1 <- temp <- data
  } else {temp.1 <- temp <- data[,-(1:(starting.col-1))]}
  
  # CP calculation
  temp.1[,] <- 0
  # 1 if lowerbound captures theta, 0 otherwise
  temp.1[,seq(1,2*len.col-1,by=2)] = (temp[,seq(1,2*len.col-1,by=2)] < theta)
  # 1 if upperbound captures theta, 0 otherwise
  temp.1[,seq(2,2*len.col,by=2)] = (temp[,seq(2,2*len.col,by=2)] >= theta)
  
  temp.2 = temp.1[,seq(1,len.col*2-1,by=2)]*temp.1[,seq(2,len.col*2,by=2)]
  names(temp.2) = CI.method
  # CP
  eval[1,] = apply(temp.2, 2,sum,na.rm=na.rm)/n.sim
  
  temp.RN = (temp.1[,seq(1,len.col*2-1,by=2)] == 1) * (temp.1[,seq(2,len.col*2,by=2)] == 0)
  temp.LN = (temp.1[,seq(1,len.col*2-1,by=2)] == 0) * (temp.1[,seq(2,len.col*2,by=2)] == 1)
  # NCP
  eval[2,] = apply(temp.LN, 2, sum, na.rm=na.rm)/n.sim
  eval[3,] = apply(temp.RN, 2, sum, na.rm=na.rm)/n.sim
  
  # NaN prob.
  eval[4,] = apply(temp.2, 2, function(x) mean(is.na(x)))
  
  # CIL]
  temp.3 = temp[,seq(2,2*len.col,by=2)] - temp[,seq(1,2*len.col-1,by=2)]
  eval[5,] = apply(temp.3, 2, mean, na.rm=na.rm)
  # ZWI
  eval[6,] = apply((temp.3 == 0), 2, sum, na.rm=na.rm)/n.sim

  
  return(round(eval,round))
  # if want to remove scientific notations,
  # apply "format(xx, scientific=F)" to each single cells, not to the whole dataframe(subject to error)
  # or manipulate the global setting: option(scipen = 999) and then return it by option(scipen = 0).
}