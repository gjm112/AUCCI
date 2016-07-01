library(ROCR)
library(mice)
library(norm)
library(ggplot2)
library(reshape2)
library(xtable)

# data step
ADdata <- read.csv("NACC_cho06172016.csv", header=T)

View(ADdata)
#### Longitudinal data to Cross-sectional (1st visit only)
# 111,223 obs -> 33,900 obs
ADdata1 <- ADdata[ADdata$NACCVNUM==1,]

### disease status: NACCETPR (1=D, 99=missing, others=ND)
### biomarkers:     NACCMMSE or CDRSUM
### covariates:     SEX, EDUC, NACCFAM, SMOKYRS, BPSYS, HRATE, PDNORMAL, 
#                   NACCGDS(depression), NACCAGEB(age at initial visit), NACCNIHR(race)
#                   NACCAGE, NACCAMD(total # of medicines), NACCBMI

### data screening
apply(ADdata1,2,unique)
ADdata1$SEX <- as.factor(ADdata1$SEX)

ADdata1$EDUC[ADdata1$EDUC==99] <- NA
ADdata1$EDUC <- as.numeric(ADdata1$EDUC)

ADdata1$NACCFAM[ADdata1$NACCFAM %in% c(9,-4)] <- NA
ADdata1$NACCFAM <- as.factor(ADdata1$NACCFAM)

ADdata1$BPSYS[ADdata1$BPSYS %in% c(888,-4)] <- NA
ADdata1$HRATE[ADdata1$HRATE %in% c(888,-4)] <- NA
ADdata1$PDNORMAL[ADdata1$PDNORMAL %in% c(8,-4)] <- NA
ADdata1$PDNORMAL <- as.factor(ADdata1$PDNORMAL)
ADdata1$NACCGDS[ADdata1$NACCGDS %in% c(88,-4)] <- NA
ADdata1$NACCMMSE[ADdata1$NACCMMSE %in% c(31:100,-4)] <- NA

# NACCETPR: disease status
ADdata1$NACCETPR[ADdata1$NACCETPR %in% c(80:100)] <- NA # missing
ADdata1$NACCETPR[ADdata1$NACCETPR %in% c(2:30)] <- 0    # symptoms other than AD
#rho (missing ratio) 48.2% missing
mean(is.na(ADdata1$NACCETPR))

# ADdata1$NACCNIHR[ADdata1$NACCNIHR == 99] <- NA
# ADdata1$NACCNIHR <- as.factor(ADdata1$NACCNIHR)
# 1=White, 2=Black or African American, 3=American Indian or Alaska Native
# 4=Native Hawaiian or Pacific Islander, 5=Asian, 6=Multiracial
# 99=Unknown or ambiguous  -> not actually missing, but unclassified!

ADdata1$NACCAMD[ADdata1$NACCAMD %in% c(-4)] <- NA # missing
ADdata1$NACCBMI[ADdata1$NACCBMI %in% c(999,-9)] <- NA # missing

apply(ADdata1,2,unique)

# glance at the data
if (FALSE) {
  tmp <- ADdata1$NACCBMI
  class(tmp); unique(tmp); levels(tmp); tmp[1:20]  
}

# reordering, dropping out unnecessary variables
ADdata1 <- ADdata1[,c("NACCID","NACCETPR","CDRSUM","CDRGLOB","SEX","NACCNIHR","SMOKYRS","NACCAGE","NACCMMSE","EDUC","NACCFAM","BPSYS","HRATE","NACCGDS","NACCAMD","NACCBMI", "PDNORMAL")]
# ID: NACCID
# disease status: NACCETPR
# biomarker candidates:CDRSUM(complete), CDRGLOB(complete), NACCMMSE(incomplete)
# complete/incomplete: NACCID~NACCAGE (complete), NACCMMSE ~ PDNORMAL (incomplete)

## randomly choosing 200 obs
set.seed(10000)
smp <- sample(1:dim(ADdata1)[1], 200)


### ADdata2 ~ ADdata5: 200 random obs, 
### 2: CDRSUM marker/ complete variables only   3: CDRSUM marker/ all variables as covariates
### 4: CDRGLOB marker/ complete variables only  5: CDRGLOB marker/ all variables as covariates

MIseed = 100

### ADdata2: with complete variables only. CDRSUM as the biomarker. randomly chosen 200 obs only
  ADdata2 <- ADdata1[smp, c("NACCID","NACCETPR","CDRSUM","SEX","NACCNIHR","SMOKYRS","NACCAGE")]
  names(ADdata2)[2] <- "diseaseR"
  names(ADdata2)[3] <- "marker"
  
  predM = 1 - diag(1, ncol(ADdata2)); predM[,1] <- 0 # ignoring the ID when predicting
  set.seed(MIseed)
  tmp <- mice(ADdata2, m=10, method="pmm", printFlag=F, predictorMatrix=predM)
  
  ADdata2.pmm <- lapply(1:10, function(x) complete(tmp, action=x))
  ADdata2.pmm.CI <- as.data.frame(t(sapply (c("Bm", "HM1", "HM2", "NW", "DL"), function(x) AUCCI.MI(ADdata2.pmm, CI.method=x)$CI)))
  apply(ADdata2.pmm.CI,1,mean)[1]   #AUC.hat
  ADdata2.pmm.CI

### ADdata3: with all variables. CDRSUM as the biomarker (and excluding CDRGLOB). randomly chosen 200 obs only
  ADdata3 <- ADdata1[smp, c("NACCID","NACCETPR","CDRSUM","SEX","NACCNIHR","SMOKYRS","NACCAGE","NACCMMSE","EDUC","NACCFAM","BPSYS","HRATE","NACCGDS","NACCAMD","NACCBMI", "PDNORMAL")]
  names(ADdata3)[2] <- "diseaseR"
  names(ADdata3)[3] <- "marker"
  
  table(ADdata3$diseaseR, useNA="always")
  mean(is.na(ADdata3$diseaseR))                     # missing coverage
  mean(as.numeric(ADdata3$diseaseR), na.rm=TRUE)  # prevalence rate of the observed

  predM = 1 - diag(1, ncol(ADdata3)); predM[,1] <- 0 # ignoring the ID when predicting
  set.seed(MIseed)
  tmp <- mice(ADdata3, m=10, method="pmm", printFlag=F, predictorMatrix=predM)
  
  ADdata3.pmm <- lapply(1:10, function(x) complete(tmp, action=x))
  ADdata3.pmm.CI <- as.data.frame(t(sapply (c("Bm", "HM1", "HM2", "NW", "DL"), function(x) AUCCI.MI(ADdata3.pmm, CI.method=x)$CI)))
  apply(ADdata3.pmm.CI,1,mean)[1]   #AUC.hat
  ADdata3.pmm.CI
  

### ADdata4: with complete variables only. CDRGLOB as the biomarker. randomly chosen 200 obs only
  ADdata4 <- ADdata1[smp, c("NACCID","NACCETPR","CDRGLOB","SEX","NACCNIHR","SMOKYRS","NACCAGE")]
  names(ADdata4)[2] <- "diseaseR"
  names(ADdata4)[3] <- "marker"
  
  predM = 1 - diag(1, ncol(ADdata4)); predM[,1] <- 0 # ignoring the ID when predicting
  set.seed(MIseed)
  tmp <- mice(ADdata4, m=10, method="pmm", printFlag=F, predictorMatrix=predM)
  
  ADdata4.pmm <- lapply(1:10, function(x) complete(tmp, action=x))
  ADdata4.pmm.CI <- as.data.frame(t(sapply (c("Bm", "HM1", "HM2", "NW", "DL"), function(x) AUCCI.MI(ADdata4.pmm, CI.method=x)$CI)))
  apply(ADdata4.pmm.CI,1,mean)[1]   #AUC.hat
  ADdata4.pmm.CI

### ADdata5: with all variables. CDRGLOB as the biomarker (and excluding CDRSUM). randomly chosen 200 obs only
  ADdata5 <- ADdata1[smp, c("NACCID","NACCETPR","CDRGLOB","SEX","NACCNIHR","SMOKYRS","NACCAGE","NACCMMSE","EDUC","NACCFAM","BPSYS","HRATE","NACCGDS","NACCAMD","NACCBMI", "PDNORMAL")]
  names(ADdata5)[2] <- "diseaseR"
  names(ADdata5)[3] <- "marker"
  
  predM = 1 - diag(1, ncol(ADdata5)); predM[,1] <- 0 # ignoring the ID when predicting
  set.seed(MIseed)
  tmp <- mice(ADdata5, m=10, method="pmm", printFlag=F, predictorMatrix=predM)
  
  ADdata5.pmm <- lapply(1:10, function(x) complete(tmp, action=x))
  ADdata5.pmm.CI <- as.data.frame(t(sapply (c("Bm", "HM1", "HM2", "NW", "DL"), function(x) AUCCI.MI(ADdata5.pmm, CI.method=x)$CI)))
  apply(ADdata5.pmm.CI,1,mean)[1]   #AUC.hat
  ADdata5.pmm.CI


### naive estimator (= 0.6435 and 0.6122)
AUC(data=ADdata2, disease="diseaseR") #CDRSUM as biomarker
AUC(data=ADdata3, disease="diseaseR") #CDRSUM as biomarker
AUC(data=ADdata4, disease="diseaseR") #CDRGLOB as biomarker
AUC(data=ADdata5, disease="diseaseR") #CDRGLOB as biomarker


### pmm, lr, norm
predM = 1 - diag(1, ncol(ADdata3)); predM[,1] <- 0 # ignoring the ID when predicting
set.seed(MIseed)
tmp <- mice(ADdata3, m=10, method="logreg", printFlag=F, predictorMatrix=predM)

set.seed(MIseed)
ADdata3.lr <- lapply(1:10, function(x) complete(tmp, action=x))
ADdata3.lr.CI <- as.data.frame(t(sapply (c("Bm", "HM1", "HM2", "NW", "DL"), function(x) AUCCI.MI(ADdata3.lr, CI.method=x)$CI)))
apply(ADdata3.lr.CI,1,mean)[1]   #AUC.hat
ADdata3.lr.CI

set.seed(MIseed)
ADdata3.NORM <- MI.norm(ADdata3[,-1], m=10, rounding="adaptive", showits=F)
ADdata3.NORM.CI <- as.data.frame(t(sapply (c("Bm", "HM1", "HM2", "NW", "DL"), function(x) AUCCI.MI(ADdata3.NORM, CI.method=x)$CI)))
apply(ADdata3.NORM.CI,1,mean)[1]   #AUC.hat
ADdata3.NORM.CI
AUCCI(data=ADdata3.NORM[[1]], disease="diseaseR", CI.method="HM1", variance=T)

ADdata3.CI <- rbind(ADdata3.pmm.CI,ADdata3.lr.CI,ADdata3.NORM.CI)

names(ADdata3.CI) = c("lb", "ub")
ADdata3.CI$CI.method = c("Bm", "HM1", "HM2", "NW", "DL")
ADdata3.CI$no = c(5,3,2,1,4)   #numbering for labels in the order of abc
ADdata3.CI$MI = rep(c("PMM","LR","NORM"), each=5)
ADdata3.CI$MI <- as.factor(ADdata3.CI$MI)
ADdata3.CI$MI <- factor(ADdata3.CI$MI,levels(ADdata3.CI$MI)[c(3,1,2)])
ADdata3.CI$length = ADdata3.CI$ub - ADdata3.CI$lb
aggregate(lb ~ MI, data=ADdata3.CI, FUN=mean)
aggregate(ub ~ MI, data=ADdata3.CI, FUN=mean)
aggregate(ub - lb ~ MI, data=ADdata3.CI, FUN=mean)
apply(ADdata3.CI[,1:2], 1, mean)

table.AD <- reshape(ADdata3.CI[,-4], idvar =c("CI.method"), timevar="MI", direction="wide")
table.AD <- rbind(table.AD, c(0,apply(table.AD[,-1],2,mean)))
table.AD[6,1] <- "average"
# table.AD[,-1] <- round(table.AD[,-1],4)
print(xtable(table.AD, digits = 4, caption="interval estimates of the AUC for the AD data"), include.rownames = FALSE)


vline.data = data.frame(AUC.hat = rep(NA,3), MI = c("PMM", "LR", "NORM"))
vline.data$AUC.hat[1] <- mean(sapply(1:10, function(x) AUC(data=ADdata3.pmm[[x]], disease="diseaseR")))
vline.data$AUC.hat[2] <- mean(sapply(1:10, function(x) AUC(data=ADdata3.lr[[x]], disease="diseaseR")))
vline.data$AUC.hat[3] <- mean(sapply(1:10, function(x) AUC(data=ADdata3.NORM[[x]], disease="diseaseR")))

AUC.naive <- AUC(data=ADdata3, disease="diseaseR")

p <- ggplot(xlim=c(0.5,1.1)) 
p +facet_grid(MI ~ .) + geom_segment(aes(x=lb, xend=ub, y=no, yend=no, color=CI.method, label=CI.method), data=ADdata3.CI) +
  geom_vline(aes(xintercept = AUC.hat),  color = "black", vline.data, show_guide = TRUE) + 
  #geom_vline(xintercept = AUC.complete,  color = "gray") +
  geom_vline(xintercept = AUC.naive, linetype="dotted", color = "blue", size=1) +
  labs(x = "confidence intervals", y = "CI methods by MI techniques") +
  theme(axis.text.y = element_blank(), axis.ticks.y =  element_blank())

ggsave("R/AD_plot.png", width = 200, height = 100, units = "mm")
