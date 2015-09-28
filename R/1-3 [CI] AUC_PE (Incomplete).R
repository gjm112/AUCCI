### Part I. CI functions
## 1-3 [CI] AUC_Point Estimate (Incomplete):  1.3.1 SeSp,  1.3.2 SeSp2AUC,   1.3.3 AUC.verif

## 1.3.1 SeSp ###########################################################################        
SeSp <- function(data, 
                 names = list("disease", "marker", c("V.1", "V.2", "V.3", "V.4", "V.5"), "R", "diseaseR"), 
                 method, cutoff=NA) {
  # data structure: D(disease), T(marker), V(covariates), R(verification), DR(diseaseR)
  # methods: full, naive, BG, MS, IPW, SP
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
  
  else if(method=="naive") {
    result = data.frame(cutoff = NA, Se = NA, Sp = NA)
    i=0; for (co in CUTOFF) {i = i+1; result[i,] = SeSp.base(T, co, D*(1-R), (1-D)*(1-R))}
  }
  
  else if(method=="BG") {
    model = as.formula(paste(names[[5]], "~", names[[2]], "*(", paste(names[[3]],collapse="+"), ")" ))
    PrD = glm(model, binomial, data=data, na.action=na.omit); PrDi = expit(predict(PrD,data))
    result = data.frame(cutoff = NA, Se = NA, Sp = NA)
    i=0; for (co in CUTOFF) {i = i+1; result[i,] = SeSp.base(T, co, PrDi, 1-PrDi)}
  }
  
  else if(method=="MS") {
    model = as.formula(paste(names[[5]], "~", names[[2]], "*(", paste(names[[3]],collapse="+"), ")" ))
    PrD = glm(model, binomial, data=data, na.action=na.omit); PrDi = expit(predict(PrD,data))
    result = data.frame(cutoff = NA, Se = NA, Sp = NA)
    i=0; for (co in CUTOFF) {i = i+1; result[i,] = SeSp.base(T, co, ((1-R)*D + R*PrDi), ((1-R)*(1-D) + R*(1-PrDi)))}
  }  
  
  else if(method=="IPW") {
    model = as.formula(paste("1-",names[[4]], "~", names[[2]], "*(", paste(names[[3]],collapse="+"), ")" ))
    PrR0 = glm(model, binomial, data=data); Pii = expit(predict(PrR0,data)) # Pii = Pr(R=0)
    result = data.frame(cutoff = NA, Se = NA, Sp = NA)
    i=0; for (co in CUTOFF) {i = i+1; result[i,] = SeSp.base(T, co, (1-R)*D*Pii^(-1), (1-R)*(1-D)*Pii^(-1))}
  } 

  else if(method=="SP") {
    model1 = as.formula(paste(names[[5]], "~", names[[2]], "*(", paste(names[[3]],collapse="+"), ")" ))
    model2 = as.formula(paste("1-",names[[4]], "~", names[[2]], "*(", paste(names[[3]],collapse="+"), ")" ))
    PrD = glm(model1, binomial, data=data, na.action=na.omit); PrDi = expit(predict(PrD, data))
    PrR0 = glm(model2, binomial, data=data); Pii = expit(predict(PrR0, data)) # Pii = Pr(R=0)
    result = data.frame(cutoff = NA, Se = NA, Sp = NA)
    i=0; for (co in CUTOFF) {i = i+1; result[i,] = SeSp.base(T, co, ((1-R)*D + (Pii-1+R)*PrDi)*Pii^(-1), ((1-R)*(1-D) + (Pii-1+R)*(1-PrDi))*Pii^(-1))}
  }
  
  return(result)
}

## 1.3.2 SeSp to AUC function ###########################################################  
SeSp2AUC = function(data) {
  Se <- c(data$Se,0)
  d.Sp <- c(data$Sp,1) - c(0, data$Sp)
  return(sum(Se*d.Sp))
}

# Example: each direct method ###########################################################
a <- SeSp(temp,,method="BG"); plot(1-a$Sp, a$Se); SeSp2AUC(a)
a <- SeSp(temp,,method="MS"); plot(1-a$Sp, a$Se); SeSp2AUC(a)
a <- SeSp(temp,,method="IPW"); plot(1-a$Sp, a$Se); SeSp2AUC(a)
a <- SeSp(temp,,method="SP"); plot(1-a$Sp, a$Se); SeSp2AUC(a)   # some perturbations happen.
a <- SeSp(temp,,method="full"); plot(1-a$Sp, a$Se); SeSp2AUC(a) 

#### 1.3.3 AUC.verif  ###################################################################  
# AUC.verif function (1.2.1 + 1.2.2 : BG, MS, IPW, SP) plus "He(direct)"
AUC.verif <- function(data,
                      names = list("disease", "marker", c("V.1", "V.2", "V.3", "V.4", "V.5"), "R", "diseaseR"), 
                      method, cutoff=NA) {
  if (method %in% c("full", "naive", "BG", "MS", "IPW", "SP")) {
    return(SeSp2AUC(SeSp(data,,method=method)))
  }
  else if (method=="He") {
    model = as.formula(paste("1-",names[[4]], "~", names[[2]], "*(", paste(names[[3]],collapse="+"), ")" ))
    PrR0 = glm(model, binomial, data=data); Pii = expit(predict(PrR0, data)) # Pii = Pr(R=0)
    marker.y = data[data$disease==1 & data$R==0, names[[2]]]
    marker.x = data[data$disease==0 & data$R==0, names[[2]]]
    Pii.y = Pii[data$disease==1 & data$R==0]
    Pii.x = Pii[data$disease==0 & data$R==0]
    n.y = length(marker.y); n.x = length(marker.x)
    A <- 0; for (i in 1:n.x) { A = A + ( (marker.y > marker.x[i]) %*% (Pii.y^(-1)) * Pii.x[i]^(-1) + ((marker.y == marker.x[i])/2) %*% (Pii.y^(-1)) * Pii.x[i]^(-1) )  }
    A = A/ sum(Pii.x^(-1)) / sum(Pii.y^(-1))
    return(A[1])
  }   
}

# Example: true value #SP is closest #####################################################
AUC(temp$marker[temp$disease==0], temp$marker[temp$disease==1])

AUC.verif(temp,,method="BG")
AUC.verif(temp,,method="MS")
AUC.verif(temp,,method="IPW")
AUC.verif(temp,,method="SP")
AUC.verif(temp,,method="full")
AUC.verif(temp,,method="naive")
AUC.verif(temp,,method="He")    # basically He is equivalent to IPW

