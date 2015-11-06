### Part I. CI functions
## 1-3 [CI] AUC_Point Estimate (Incomplete):  1.3.1 SeSp,  1.3.2 SeSp2AUC,   1.3.3 AUC.verif

## 1.3.1 SeSp ###########################################################################        
SeSp <- function(data, 
                 names = list("disease", "marker", c("V.1", "V.2", "V.3", "V.4", "V.5"), "R", "diseaseR"), 
                 CI.method, cutoff=NA) {
  # data structure: D(disease), T(marker), V(covariates), R(verification), DR(diseaseR)
  # CI.methods: full, naive, BG, MS, IPW, SP
  # single.cutoff
  # cutoff.approx(number of cutoff): when n(T) is too large, set up the number of cutoffs to approximate
  D = data[,names[[1]]]
  T = data[,names[[2]]]
  V = data[,names[[3]]]
  R = data[,names[[4]]]
  DR = data[,names[[5]]]
  if (is.na(cutoff)) {CUTOFF=sort(T)} else {CUTOFF=cutoff}
  
  SeSp.base = function(T, cutoff, weight.Se, weight.Sp) {
    Se = sum((T>=cutoff)*weight.Se)/sum(weight.Se)
    Sp = sum((T<cutoff)*weight.Sp)/sum(weight.Sp)
    return(c(cutoff, Se, Sp))
  }
  
  result = list()
  
  if("full" %in% CI.method) {
    result[["full"]] = data.frame(cutoff = NA, Se = NA, Sp = NA)
    i=0; for (co in CUTOFF) {i = i+1; result[["full"]][i,] = SeSp.base(T, co, D, 1-D)}
  }
  
  if("naive" %in% CI.method) {
    result[["naive"]] = data.frame(cutoff = NA, Se = NA, Sp = NA)
    i=0; for (co in CUTOFF) {i = i+1; result[["naive"]][i,] = SeSp.base(T, co, D*(1-R), (1-D)*(1-R))}
  }
  
  if("BG" %in% CI.method) {
    model = as.formula(paste(names[[5]], "~", names[[2]], "*(", paste(names[[3]],collapse="+"), ")" ))
    PrD = glm(model, binomial, data=data, na.action=na.omit); PrDi = expit(predict(PrD,data))
    result[["BG"]] = data.frame(cutoff = NA, Se = NA, Sp = NA)
    i=0; for (co in CUTOFF) {i = i+1; result[["BG"]][i,] = SeSp.base(T, co, PrDi, 1-PrDi)}
  }
  
  if("MS" %in% CI.method) {
    model = as.formula(paste(names[[5]], "~", names[[2]], "*(", paste(names[[3]],collapse="+"), ")" ))
    PrD = glm(model, binomial, data=data, na.action=na.omit); PrDi = expit(predict(PrD,data))
    result[["MS"]] = data.frame(cutoff = NA, Se = NA, Sp = NA)
    i=0; for (co in CUTOFF) {i = i+1; result[["MS"]][i,] = SeSp.base(T, co, ((1-R)*D + R*PrDi), ((1-R)*(1-D) + R*(1-PrDi)))}
  }  
  
  if("IPW" %in% CI.method) {
    model = as.formula(paste("1-",names[[4]], "~", names[[2]], "*(", paste(names[[3]],collapse="+"), ")" ))
    PrR0 = glm(model, binomial, data=data); Pii = expit(predict(PrR0,data)) # Pii = Pr(R=0)
    result[["IPW"]] = data.frame(cutoff = NA, Se = NA, Sp = NA)
    i=0; for (co in CUTOFF) {i = i+1; result[["IPW"]][i,] = SeSp.base(T, co, (1-R)*D*Pii^(-1), (1-R)*(1-D)*Pii^(-1))}
  } 

  if("SP" %in% CI.method) {
    model1 = as.formula(paste(names[[5]], "~", names[[2]], "*(", paste(names[[3]],collapse="+"), ")" ))
    model2 = as.formula(paste("1-",names[[4]], "~", names[[2]], "*(", paste(names[[3]],collapse="+"), ")" ))
    PrD = glm(model1, binomial, data=data, na.action=na.omit); PrDi = expit(predict(PrD, data))
    PrR0 = glm(model2, binomial, data=data); Pii = expit(predict(PrR0, data)) # Pii = Pr(R=0)
    result[["SP"]] = data.frame(cutoff = NA, Se = NA, Sp = NA)
    i=0; for (co in CUTOFF) {i = i+1; result[["SP"]][i,] = SeSp.base(T, co, ((1-R)*D + (Pii-1+R)*PrDi)*Pii^(-1), ((1-R)*(1-D) + (Pii-1+R)*(1-PrDi))*Pii^(-1))}
  }
  
  return(result)
}

## 1.3.2 SeSp to AUC function ###########################################################  
SeSp2AUC = function(data) {
  if (class(data) != "list") {stop("Data is not a list")}
  len = length(data)
  result = data.frame(t(rep(NA,len)))
  colnames(result) <- names(data)
  for (i in 1:len){
    Se <- c(data[[i]]$Se,0)
    d.Sp <- c(data[[i]]$Sp,1) - c(0, data[[i]]$Sp)
    result[,i] = sum(Se*d.Sp)
  }
  return(result)
}

#### 1.3.3 AUC.verif  ###################################################################  
# AUC.verif function (1.2.1 + 1.2.2 : BG, MS, IPW, SP) plus "He(direct)"
AUC.verif <- function(data,
                      names = list("disease", "marker", c("V.1", "V.2", "V.3", "V.4", "V.5"), "R", "diseaseR"), 
                      CI.method, cutoff=NA) {
  len = length(CI.method)
  CI.SeSp = CI.method[CI.method %in% c("full", "naive", "BG", "MS", "IPW", "SP")]
  if (length(CI.SeSp)>0) {result = SeSp2AUC(SeSp(data,,CI.method=CI.SeSp))} else {result = data.frame(He=NA)}
  if ("He" %in% CI.method) {
    model = as.formula(paste("1-",names[[4]], "~", names[[2]], "*(", paste(names[[3]],collapse="+"), ")" ))
    PrR0 = glm(model, binomial, data=data); Pii = expit(predict(PrR0, data)) # Pii = Pr(R=0)
    marker.y = data[data$disease==1 & data$R==0, names[[2]]]
    marker.x = data[data$disease==0 & data$R==0, names[[2]]]
    Pii.y = Pii[data$disease==1 & data$R==0]
    Pii.x = Pii[data$disease==0 & data$R==0]
    n.y = length(marker.y); n.x = length(marker.x)
    if (n.x * n.y == 0) { A = NA } else {         # in case of no x or no y in resampling, error occurs. NAing it.
      A <- 0; for (i in 1:n.x) { A = A + ( (marker.y > marker.x[i]) %*% (Pii.y^(-1)) * Pii.x[i]^(-1) + ((marker.y == marker.x[i])/2) %*% (Pii.y^(-1)) * Pii.x[i]^(-1) )  }
      A = A/ sum(Pii.x^(-1)) / sum(Pii.y^(-1))      
    }
    result$He = A[1]
  }
  return(result)
}