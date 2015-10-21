### Part II. Simulation
## 2-2 [SIM] Evaluator: 2.2.1 CI.evaluator ################################################

## 2.2.1 CI.evaluator #####################################################################
CI.evaluator = function(data, param = 2, est = 3, n.i,
                        CI.method = sub(".lb","",matrix(names(data[[est]])[-1],nrow=2)[1,]),
                        na.rm = TRUE) {
  msmt = c("CP","LNCP","RNCP","CIL", "ZWI")
  
  # data formatting
  len.col = length(CI.method); len.row = length(msmt)
  eval = as.data.frame(matrix(NA,len.row,len.col))
  names(eval) = CI.method
  rownames(eval) = msmt
  
  # informations
  theta = data[[param]]$theta
  
  # CP calculation
  temp.1 <- temp <- data[[est]][[n.i]][,-1];
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
  
  return(eval)
}