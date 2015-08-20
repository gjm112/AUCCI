### AUCCI Part I. CI functions
## 1-2. CI functions with verification bias

####### modify SeSp as a function of cutoff "vector"!!!
####### what is the estimate of pi? in IPW and SP methods

cutoff=seq(min(data[,names[[2]]]),max(data[,names[[2]]]),length.out=1000

           
## 1.2.1 SeSp function           
SeSp <- function(data, 
                 names = list("disease", "marker", c("V.1", "V.2", "V.3", "V.4", "V.5"), "R", "diseaseR"), 
                 method, cutoff=NA) {
  # data structure: D(disease), T(marker), V(covariates), R(verification), DR(diseaseR)
  # methods: BG, MS, IPW, SP
  # single.cutoff
  # cutoff.approx(number of cutoff): when n(T) is too large, set up the number of cutoffs to approximate
  D = data[,names[[1]]]
  T = data[,names[[2]]]
  V = data[,names[[3]]]
  R = data[,names[[4]]]
  DR = data[,names[[5]]]
  if (is.na(cutoff)) {CUTOFF=sort(T)} else {CUTOFF=cutoff}

  
  expit = function(x) { exp(x) / (exp(x)+1) }
  SeSp.base = function(T, cutoff, weight.Se, weight.Sp) {
    Se = sum((T>=cutoff)*weight.Se)/sum(weight.Se)
    Sp = sum((T<cutoff)*weight.Sp)/sum(weight.Sp)
    return(c(cutoff, Se, Sp))
  }
  
  if(method=="full") {
    result = data.frame(cutoff = NA, Se = NA, Sp = NA)
    i=0; for (co in CUTOFF) {i = i+1; result[i,] = SeSp.base(T, co, D, 1-D)}
  }
  
  if(method=="naive") {
    result = data.frame(cutoff = NA, Se = NA, Sp = NA)
    i=0; for (co in CUTOFF) {i = i+1; result[i,] = SeSp.base(T, co, D*(1-R), (1-D)*(1-R))}
  }
  
  if(method=="BG") {
    model = as.formula(paste(names[[5]], "~", names[[2]], "*(", paste(names[[3]],collapse="+"), ")" ))
    PrD = glm(model, binomial, data=data, na.action=na.omit); PrDi = expit(predict(PrD,data))
    result = data.frame(cutoff = NA, Se = NA, Sp = NA)
    i=0; for (co in CUTOFF) {i = i+1; result[i,] = SeSp.base(T, co, PrDi, 1-PrDi)}
  }
  
  if(method=="MS") {
    model = as.formula(paste(names[[5]], "~", names[[2]], "*(", paste(names[[3]],collapse="+"), ")" ))
    PrD = glm(model, binomial, data=data, na.action=na.omit); PrDi = expit(predict(PrD,data))
    result = data.frame(cutoff = NA, Se = NA, Sp = NA)
    i=0; for (co in CUTOFF) {i = i+1; result[i,] = SeSp.base(T, co, ((1-R)*D + R*PrDi), ((1-R)*(1-D) + R*(1-PrDi)))}
  }  
  
  if(method=="IPW") {
    model = as.formula(paste("1-",names[[4]], "~", names[[2]], "*(", paste(names[[3]],collapse="+"), ")" ))
    PrR0 = glm(model, binomial, data=data); Pii = expit(predict(PrR0,data)) # Pii = Pr(R=0)
    result = data.frame(cutoff = NA, Se = NA, Sp = NA)
    i=0; for (co in CUTOFF) {i = i+1; result[i,] = SeSp.base(T, co, (1-R)*D*Pii^(-1), (1-R)*(1-D)*Pii^(-1))}
  } 

  if(method=="SP") {
    model1 = as.formula(paste(names[[5]], "~", names[[2]], "*(", paste(names[[3]],collapse="+"), ")" ))
    model2 = as.formula(paste("1-",names[[4]], "~", names[[2]], "*(", paste(names[[3]],collapse="+"), ")" ))
    PrD = glm(model1, binomial, data=data, na.action=na.omit); PrDi = expit(predict(PrD, data))
    PrR0 = glm(model2, binomial, data=data); Pii = expit(predict(PrR0, data)) # Pii = Pr(R=0)
    result = data.frame(cutoff = NA, Se = NA, Sp = NA)
    i=0; for (co in CUTOFF) {i = i+1; result[i,] = SeSp.base(T, co, ((1-R)*D + (Pii-1+R)*PrDi)*Pii^(-1), ((1-R)*(1-D) + (Pii-1+R)*(1-PrDi))*Pii^(-1))}
  }  
  
  return(result)
}

#### 1.2.2 SeSp to AUC function
SeSp2AUC = function(data) {
#  Se.ave <- (c(1,data$Se) + c(data$Se,0))/2
  Se.ave <- c(data$Se,0)
  Sp.dif <- c(data$Sp,1) - c(0, data$Sp)
  return(sum(Se.ave*Sp.dif))
}

# each direct method
a <- SeSp(temp,,method="BG"); plot(1-a$Sp, a$Se); SeSp2AUC(a)
a <- SeSp(temp,,method="MS"); plot(1-a$Sp, a$Se); SeSp2AUC(a)
a <- SeSp(temp,,method="IPW"); plot(1-a$Sp, a$Se); SeSp2AUC(a)
a <- SeSp(temp,,method="SP"); plot(1-a$Sp, a$Se); SeSp2AUC(a)   # some perturbations happen.
a <- SeSp(temp,,method="full"); plot(1-a$Sp, a$Se); SeSp2AUC(a) 
# true value #SP is closest
AUC(temp$marker[temp$disease==0], temp$marker[temp$disease==1])
