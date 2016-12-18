library(ROCR)
library(mice)
library(norm)
library(ggplot2)
library(reshape2)
library(xtable)

# data step
ADdata <- read.csv("NACC_cho08192016.csv", header=T)

#### Longitudinal data to Cross-sectional (1st visit only)
# 111,223 obs -> 33,900 obs
ADdata1 <- ADdata[ADdata$NACCVNUM==1,]

### disease status: NACCETPR (1=D, 99=missing, others=ND)
### biomarkers:     NACCMMSE or CDRSUM
### covariates:     SEX, EDUC, NACCFAM, SMOKYRS, BPSYS, HRATE, PDNORMAL, 
#                   NACCGDS(depression), NACCAGEB(age at initial visit), NACCNIHR(race)
#                   NACCAGE, NACCAMD(total # of medicines), NACCBMI

### data cleansing
ADdata1$SEX <- as.factor(ADdata1$SEX)

ADdata1$EDUC[ADdata1$EDUC==99] <- NA
ADdata1$EDUC <- as.numeric(ADdata1$EDUC)

ADdata1$NACCFAM[ADdata1$NACCFAM %in% c(9,-4)] <- NA
ADdata1$NACCFAM <- as.factor(ADdata1$NACCFAM)

ADdata1$BPSYS[ADdata1$BPSYS %in% c(888,-4)] <- NA
ADdata1$HRATE[ADdata1$HRATE %in% c(888,-4)] <- NA

ADdata1$PDNORMAL[ADdata1$PDNORMAL %in% c(8,-4)] <- NA
ADdata1$PDNORMAL <- as.factor(ADdata1$PDNORMAL)

ADdata1$SMOKYRS[ADdata1$SMOKYRS %in% c(88,99,-4)] <- NA

ADdata1$NACCNIHR[ADdata1$NACCNIHR %in% c(88,99,-4)] <- NA
# Forcing multicategorical variable into binary
ADdata1$NACCNIHR2 <- 0
ADdata1$NACCNIHR2[ADdata1$NACCNIHR==1] <- 1
ADdata1$NACCNIHR2[is.na(ADdata1$NACCNIHR)] <- NA
ADdata1$NACCNIHR2 <- as.factor(ADdata1$NACCNIHR2)
ADdata1$NACCNIHR <- NULL

ADdata1$NACCGDS[ADdata1$NACCGDS %in% c(88,-4)] <- NA
ADdata1$NACCMMSE[ADdata1$NACCMMSE %in% c(31:100,-4)] <- NA

# NACCETPR: supplementary disease info
ADdata1$NACCETPR[ADdata1$NACCETPR %in% c(80:100)] <- NA # missing
ADdata1$NACCETPR[ADdata1$NACCETPR %in% c(2:30)] <- 0    # symptoms other than AD
ADdata1$NACCETPR <- as.factor(ADdata1$NACCETPR)

# NPADNC: disease status - golden standard test
ADdata1$NPADNC[ADdata1$NPADNC %in% c(-4,8,9)] <- NA # missing
ADdata1$disease <- ADdata1$NPADNC
ADdata1$disease[ADdata1$disease %in% c(1,2,3)] <- 1   # AD
ADdata1$disease <- as.factor(ADdata1$disease)
#rho (missing ratio) 97.1% missing
mean(is.na(ADdata1$disease))
# D 844 (2.5%), ND 135 (0.4%), missing 32921 (97%)
table(ADdata1$disease,useNA="always")
# P(D|obs)=86%
mean(as.numeric(ADdata1$disease)-1,na.rm=TRUE)


# ADdata1$NACCNIHR[ADdata1$NACCNIHR == 99] <- NA
# ADdata1$NACCNIHR <- as.factor(ADdata1$NACCNIHR)
# 1=White, 2=Black or African American, 3=American Indian or Alaska Native
# 4=Native Hawaiian or Pacific Islander, 5=Asian, 6=Multiracial
# 99=Unknown or ambiguous  -> not actually missing, but unclassified!

ADdata1$NACCAMD[ADdata1$NACCAMD %in% c(-4)] <- NA # missing
ADdata1$NACCBMI[ADdata1$NACCBMI %in% c(999,888,-9,-4)] <- NA # missing
ADdata1$NACCBMI[ADdata1$NACCBMI >=888] <- NA # missing (to handle 888.89999)
# reordering, dropping out unnecessary variables
# ADdata1 <- ADdata1[,c("NACCID","NACCETPR","CDRSUM","CDRGLOB","SEX","NACCNIHR2","SMOKYRS","NACCAGE","NACCMMSE","EDUC","NACCFAM","BPSYS","HRATE","NACCGDS","NACCAMD","NACCBMI", "PDNORMAL")]
# ID: NACCID
# disease status: NACCETPR
# biomarker candidates:CDRSUM(complete), CDRGLOB(complete), NACCMMSE(incomplete)
# complete/incomplete: NACCID~NACCAGE (complete), NACCMMSE ~ PDNORMAL (incomplete)

## randomly choosing one of the centers (NACCADC)
table(ADdata1$NACCADC, ADdata1$disease,useNA="always")

# an ADC center example (4347)
table(ADdata1$disease[ADdata1$NACCADC==4347],useNA="always")
### ADdata2: ADC==4347, CDRSUM marker/ all variables as covariates
ADdata2 <- ADdata1[ADdata1$NACCADC==4347, c("NACCID","disease","CDRSUM","SEX","NACCETPR", "NACCNIHR2","SMOKYRS","NACCAGE","NACCMMSE","EDUC","NACCFAM","BPSYS","HRATE","NACCGDS","NACCAMD","NACCBMI", "PDNORMAL")]
#table(ADdata2$disease, useNA="always")

names(ADdata2)[2] <- "diseaseR"
names(ADdata2)[3] <- "marker"

table(ADdata2$diseaseR, useNA="always")
mean(is.na(ADdata2$diseaseR))                     # missing coverage 93.4%
mean(as.numeric(ADdata2$diseaseR), na.rm=TRUE)-1  # prevalence rate of the observed 83.6%

## 0. naive estimator (= 0.589)
AUC(data=ADdata2, disease="diseaseR") #CDRSUM as biomarker


## 1. PMM
predM = 1 - diag(1, ncol(ADdata2)); predM[,1] <- 0 # ignoring the ID when predicting
set.seed(MIseed)
tmp <- mice(ADdata2, m=10, method="pmm", printFlag=F, predictorMatrix=predM)

ADdata2.pmm <- lapply(1:10, function(x) complete(tmp, action=x))
ADdata2.pmm.CI <- as.data.frame(t(sapply (c("Bm", "HM1", "HM2", "NW", "DL"), function(x) AUCCI.MI(ADdata2.pmm, CI.method=x)$CI)))
apply(ADdata2.pmm.CI,1,mean)[1]   #AUC.hat
ADdata2.pmm.CI


## 2. LR
predM = 1 - diag(1, ncol(ADdata2)); predM[,1] <- 0 # ignoring the ID when predicting
set.seed(MIseed)
method.tmp = rep("pmm",17); method.tmp[c(2,4,5,11,17)] <- "logreg"
tmp <- mice(ADdata2, m=10, method=method.tmp, printFlag=F, predictorMatrix=predM)
ADdata2.lr <- lapply(1:10, function(x) complete(tmp, action=x))
ADdata2.lr.CI <- as.data.frame(t(sapply (c("Bm", "HM1", "HM2", "NW", "DL"), function(x) AUCCI.MI(ADdata2.lr, CI.method=x)$CI)))
apply(ADdata2.lr.CI,1,mean)[1]   #AUC.hat
ADdata2.lr.CI


## 3. NORM
#As numeric for NORM
ADdata2$diseaseR <- as.numeric(as.character(ADdata2$diseaseR))
ADdata2$SEX <- as.numeric(as.character(ADdata2$SEX))
ADdata2$NACCFAM <- as.numeric(as.character(ADdata2$NACCFAM))
ADdata2$PDNORMAL <- as.numeric(as.character(ADdata2$PDNORMAL))
ADdata2$NACCNIHR2 <- as.numeric(as.character(ADdata2$NACCNIHR2))

set.seed(MIseed)
ADdata2.NORM <- MI.norm2(data=ADdata2[,-1], m=10, rounding="adaptive", showits=F) #NACCNIHR2(17th col) instead of NACCNIHR(5th)
ADdata2.NORM.CI <- as.data.frame(t(sapply (c("Bm", "HM1", "HM2", "NW", "DL"), function(x) AUCCI.MI(ADdata2.NORM, CI.method=x)$CI)))
apply(ADdata2.NORM.CI,1,mean)[1]   #AUC.hat
ADdata2.NORM.CI
#AUCCI(data=ADdata2.NORM[[1]], disease="diseaseR", CI.method="HM1", variance=T)


## Combine all MIs
ADdata2.CI <- rbind(ADdata2.pmm.CI,ADdata2.lr.CI,ADdata2.NORM.CI)
ADdata2.CI

names(ADdata2.CI) = c("lb", "ub")
ADdata2.CI$CI.method = c("Bm", "HM1", "HM2", "NW", "DL")
ADdata2.CI$no = c(5,3,2,1,4)   #numbering for labels in the order of abc
ADdata2.CI$MI = rep(c("PMM","LR","NORM"), each=5)
ADdata2.CI$MI <- as.factor(ADdata2.CI$MI)
ADdata2.CI$MI <- factor(ADdata2.CI$MI,levels(ADdata2.CI$MI)[c(3,1,2)])
ADdata2.CI$length = ADdata2.CI$ub - ADdata2.CI$lb

## descriptive stats
summary(ADdata2)


## table for LaTex
table.AD <- reshape(ADdata2.CI[,-4], idvar =c("CI.method"), timevar="MI", direction="wide")
table.AD <- rbind(table.AD, c(0,apply(table.AD[,-1],2,mean)))
table.AD[6,1] <- "average"
print(xtable(table.AD, digits = 4, caption="interval estimates of the AUC for the AD data"), include.rownames = FALSE)


## plots
vline.data = data.frame(AUC.hat = rep(NA,3), MI = c("PMM", "LR", "NORM"))
vline.data$AUC.hat[1] <- mean(sapply(1:10, function(x) AUC(data=ADdata2.pmm[[x]], disease="diseaseR")))
vline.data$AUC.hat[2] <- mean(sapply(1:10, function(x) AUC(data=ADdata2.lr[[x]], disease="diseaseR")))
vline.data$AUC.hat[3] <- mean(sapply(1:10, function(x) AUC(data=ADdata2.NORM[[x]], disease="diseaseR")))

AUC.naive <- AUC(data=ADdata2, disease="diseaseR")

p <- ggplot(xlim=c(min(ADdata2.CI$lb)-.05,0.9))
p + facet_grid(MI ~ .) + 
  geom_segment(aes(x=lb, xend=ub, y=no, yend=no, color="grey", label=NULL), size=1.5,  
               data=ADdata2.CI) +
  geom_point(aes(AUC.hat, 1), color = "black", size = 2, data=vline.data) +
  geom_point(aes(AUC.hat, 2), color = "black", size = 2, data=vline.data) +
  geom_point(aes(AUC.hat, 3), color = "black", size = 2, data=vline.data) +
  geom_point(aes(AUC.hat, 4), color = "black", size = 2, data=vline.data) +
  geom_point(aes(AUC.hat, 5), color = "black", size = 2, data=vline.data) +
  #geom_vline(aes(xintercept = AUC.hat), color = "black", vline.data) + 
  geom_vline(xintercept = AUC.naive, linetype="dotted", color = "blue", size=1) +
  labs(x = "confidence intervals", y = "CI methods by MI techniques") +
  theme(axis.text.y = element_blank(), axis.ticks.y =  element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "grey80"),
        legend.position="none")+
  #scale_color_discrete(name="CI methods") +
  geom_text(data=ADdata2.CI, x=0.77, y=ADdata2.CI$no, aes(label=CI.method), show.legend = FALSE)+
  xlim(c(min(ADdata2.CI$lb)-.05,0.9))


ggsave("R/plot_AD.png", width = 200, height = 100, units = "mm")

