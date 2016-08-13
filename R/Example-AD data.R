library(ROCR)
library(mice)
library(norm)
library(ggplot2)
library(reshape2)
library(xtable)

# data step
ADdata <- read.csv("NACC_cho06172016.csv", header=T)

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

# NACCETPR: disease status
ADdata1$NACCETPR[ADdata1$NACCETPR %in% c(80:100)] <- NA # missing
ADdata1$NACCETPR[ADdata1$NACCETPR %in% c(2:30)] <- 0    # symptoms other than AD
ADdata1$NACCETPR <- as.factor(ADdata1$NACCETPR)
#rho (missing ratio) 48.2% missing
mean(is.na(ADdata1$NACCETPR))

# ADdata1$NACCNIHR[ADdata1$NACCNIHR == 99] <- NA
# ADdata1$NACCNIHR <- as.factor(ADdata1$NACCNIHR)
# 1=White, 2=Black or African American, 3=American Indian or Alaska Native
# 4=Native Hawaiian or Pacific Islander, 5=Asian, 6=Multiracial
# 99=Unknown or ambiguous  -> not actually missing, but unclassified!

ADdata1$NACCAMD[ADdata1$NACCAMD %in% c(-4)] <- NA # missing
ADdata1$NACCBMI[ADdata1$NACCBMI %in% c(999,-9)] <- NA # missing

# reordering, dropping out unnecessary variables
ADdata1 <- ADdata1[,c("NACCID","NACCETPR","CDRSUM","CDRGLOB","SEX","NACCNIHR2","SMOKYRS","NACCAGE","NACCMMSE","EDUC","NACCFAM","BPSYS","HRATE","NACCGDS","NACCAMD","NACCBMI", "PDNORMAL")]
# ID: NACCID
# disease status: NACCETPR
# biomarker candidates:CDRSUM(complete), CDRGLOB(complete), NACCMMSE(incomplete)
# complete/incomplete: NACCID~NACCAGE (complete), NACCMMSE ~ PDNORMAL (incomplete)

## randomly choosing 200 obs
set.seed(10000)
smp <- sample(1:dim(ADdata1)[1], 200)


### ADdata1: 200 random obs, CDRSUM marker/ all variables as covariates
ADdata1 <- ADdata1[smp, c("NACCID","NACCETPR","CDRSUM","SEX","NACCNIHR2","SMOKYRS","NACCAGE","NACCMMSE","EDUC","NACCFAM","BPSYS","HRATE","NACCGDS","NACCAMD","NACCBMI", "PDNORMAL")]
names(ADdata1)[2] <- "diseaseR"
names(ADdata1)[3] <- "marker"

table(ADdata1$diseaseR, useNA="always")
mean(is.na(ADdata1$diseaseR))                     # missing coverage
mean(as.numeric(ADdata1$diseaseR), na.rm=TRUE)-1  # prevalence rate of the observed

## 0. naive estimator (= 0.67239)
AUC(data=ADdata1, disease="diseaseR") #CDRSUM as biomarker


## 1. PMM
predM = 1 - diag(1, ncol(ADdata1)); predM[,1] <- 0 # ignoring the ID when predicting
set.seed(MIseed)
tmp <- mice(ADdata1, m=10, method="pmm", printFlag=F, predictorMatrix=predM)

ADdata1.pmm <- lapply(1:10, function(x) complete(tmp, action=x))
ADdata1.pmm.CI <- as.data.frame(t(sapply (c("Bm", "HM1", "HM2", "NW", "DL"), function(x) AUCCI.MI(ADdata1.pmm, CI.method=x)$CI)))
apply(ADdata1.pmm.CI,1,mean)[1]   #AUC.hat
ADdata1.pmm.CI


## 2. LR
predM = 1 - diag(1, ncol(ADdata1)); predM[,1] <- 0 # ignoring the ID when predicting
set.seed(MIseed)
method.tmp = rep("pmm",16); method.tmp[c(2,4,5,10,16)] <- "logreg"
tmp <- mice(ADdata1, m=10, method=method.tmp, printFlag=F, predictorMatrix=predM)
ADdata1.lr <- lapply(1:10, function(x) complete(tmp, action=x))
ADdata1.lr.CI <- as.data.frame(t(sapply (c("Bm", "HM1", "HM2", "NW", "DL"), function(x) AUCCI.MI(ADdata1.lr, CI.method=x)$CI)))
apply(ADdata1.lr.CI,1,mean)[1]   #AUC.hat
ADdata1.lr.CI


## 3. NORM
#As numeric for NORM
ADdata1$diseaseR <- as.numeric(as.character(ADdata1$diseaseR))
ADdata1$SEX <- as.numeric(as.character(ADdata1$SEX))
ADdata1$NACCFAM <- as.numeric(as.character(ADdata1$NACCFAM))
ADdata1$PDNORMAL <- as.numeric(as.character(ADdata1$PDNORMAL))
ADdata1$NACCNIHR2 <- as.numeric(as.character(ADdata1$NACCNIHR2))

set.seed(MIseed)
ADdata1.NORM <- MI.norm2(data=ADdata1[,-1], m=10, rounding="adaptive", showits=F) #NACCNIHR2(17th col) instead of NACCNIHR(5th)
ADdata1.NORM.CI <- as.data.frame(t(sapply (c("Bm", "HM1", "HM2", "NW", "DL"), function(x) AUCCI.MI(ADdata1.NORM, CI.method=x)$CI)))
apply(ADdata1.NORM.CI,1,mean)[1]   #AUC.hat
ADdata1.NORM.CI
#AUCCI(data=ADdata1.NORM[[1]], disease="diseaseR", CI.method="HM1", variance=T)


## Combine all MIs
ADdata1.CI <- rbind(ADdata1.pmm.CI,ADdata1.lr.CI,ADdata1.NORM.CI)

names(ADdata1.CI) = c("lb", "ub")
ADdata1.CI$CI.method = c("Bm", "HM1", "HM2", "NW", "DL")
ADdata1.CI$no = c(5,3,2,1,4)   #numbering for labels in the order of abc
ADdata1.CI$MI = rep(c("PMM","LR","NORM"), each=5)
ADdata1.CI$MI <- as.factor(ADdata1.CI$MI)
ADdata1.CI$MI <- factor(ADdata1.CI$MI,levels(ADdata1.CI$MI)[c(3,1,2)])
ADdata1.CI$length = ADdata1.CI$ub - ADdata1.CI$lb


## table for LaTex
table.AD <- reshape(ADdata1.CI[,-4], idvar =c("CI.method"), timevar="MI", direction="wide")
table.AD <- rbind(table.AD, c(0,apply(table.AD[,-1],2,mean)))
table.AD[6,1] <- "average"
print(xtable(table.AD, digits = 4, caption="interval estimates of the AUC for the AD data"), include.rownames = FALSE)


## plots
vline.data = data.frame(AUC.hat = rep(NA,3), MI = c("PMM", "LR", "NORM"))
vline.data$AUC.hat[1] <- mean(sapply(1:10, function(x) AUC(data=ADdata3.pmm[[x]], disease="diseaseR")))
vline.data$AUC.hat[2] <- mean(sapply(1:10, function(x) AUC(data=ADdata3.lr[[x]], disease="diseaseR")))
vline.data$AUC.hat[3] <- mean(sapply(1:10, function(x) AUC(data=ADdata3.NORM[[x]], disease="diseaseR")))

AUC.naive <- AUC(data=ADdata1, disease="diseaseR")

p <- ggplot(xlim=c(0.4,1.0))
p + facet_grid(MI ~ .) + 
  geom_segment(aes(x=lb, xend=ub, y=no, yend=no, color="grey", label=NULL), size=1.5,  
               data=ADdata1.CI) +
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
  geom_text(data=ADdata1.CI, x=0.77, y=ADdata1.CI$no, aes(label=CI.method), show.legend = FALSE) +
  xlim(c(0.4,1.0))


ggsave("R/plot_AD.png", width = 200, height = 100, units = "mm")

