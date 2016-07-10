## library
library(reshape); library(reshape2)
library(xtable)
MAE = function(CP.dev , theta=.95) {mean(abs(CP.dev - theta))}
MSE = function(CP.dev) {sqrt(sum((CP.dev-theta)^2))/(length(CP.dev)-1)}
Wald <- levels(agreg.table.dev$CI.method)[c(1:3,5,8)]
MIset <- c("complete","pmm","logreg", "adaptive")
measurement <- c("CP", "LNCP", "RNCP", "CIL", "NaN", "ZWI")
measurement2 <- c(measurement,"CP.MAE")[c(1,7,2,3,4)]

## combining evaluations in a sheet
msmt = c("CP","LNCP","RNCP","NaN","CIL", "ZWI")
mi.names=c("complete","pmm","logreg","simple","coinflip","adaptive")
n = c(200,100,50)
d1 = 1:length(param1$alphabet$theta); d2 = 1:dim(param1$gamma)[1]; d3 = 1:length(n); d4 = 1:length(mi.names)


# data formatting
len.col = length(CI.methods) + 6   # 5 for th, ph, rh, n, and mi
len.row = length(msmt)*max(d1)*max(d2)*max(d3)*max(d4)
result = as.data.frame(matrix(NA,len.row,len.col))
names(result) = c("theta","phi","rho","n","MI","measure",CI.methods)
rownames(result) = NULL

index=1; increment=length(msmt)
pointest = as.data.frame(matrix(NA,288,6)) # 288 = theta(4)*phi(2)*rho(2)*n(3)*MI(6), 6: th,ph,rh,n,MI,AUC.hat
names(pointest) = c("AUC","phi","rho","n","MI","AUC.hat")
mmmm = 1 # row index for pointest

for (i in d1) {         #i: theta*phi
  th = round(param1$alphabet$theta[i],4); ph = param1$alphabet$phi[i]
  for (j in d2) {       #j: rho
    rh = param1$gamma$rho[j]
    tmp <- readRDS(paste0("R/Simdata2/sim_data-Jul01-",i,j,".rds"))
    # CI lengths to be restricted within [0,1]
    for (h in d3){
      # truncation
      tmp[[h]]$est.com[tmp[[h]]$est.com>1] = 1 ; tmp[[h]]$est.com[tmp[[h]]$est.com<0] = 0
      tmp[[h]]$est.MI$pmm[tmp[[h]]$est.MI$pmm>1] = 1 ; tmp[[h]]$est.MI$pmm[tmp[[h]]$est.MI$pmm<0] = 0
      tmp[[h]]$est.MI$logreg[tmp[[h]]$est.MI$logreg>1] = 1 ; tmp[[h]]$est.MI$logreg[tmp[[h]]$est.MI$logreg<0] = 0
      tmp[[h]]$est.MI$simple[tmp[[h]]$est.MI$simple>1] = 1 ; tmp[[h]]$est.MI$simple[tmp[[h]]$est.MI$simple<0] = 0
      tmp[[h]]$est.MI$coinflip[tmp[[h]]$est.MI$coinflip>1] = 1 ; tmp[[h]]$est.MI$coinflip[tmp[[h]]$est.MI$coinflip<0] = 0
      tmp[[h]]$est.MI$adaptive[tmp[[h]]$est.MI$adaptive>1] = 1 ; tmp[[h]]$est.MI$adaptive[tmp[[h]]$est.MI$adaptive<0] = 0
      # getting average and then rounding to 4 decimal points
      lgth = length(CI.methods); lbcol = seq(2,2*lgth+1,2); ubcol = seq(3,2*lgth+1,2)
      tmp[[h]]$eval[[1]][5,] = round(apply(tmp[[h]]$est.com[,ubcol] - tmp[[h]]$est.com[,lbcol], 2, mean, na.rm=TRUE),4)
      tmp[[h]]$eval[[2]][5,] = round(apply(tmp[[h]]$est.MI$pmm[,ubcol] - tmp[[h]]$est.MI$pmm[,lbcol], 2, mean, na.rm=TRUE),4)
      tmp[[h]]$eval[[3]][5,] = round(apply(tmp[[h]]$est.MI$logreg[,ubcol] - tmp[[h]]$est.MI$logreg[,lbcol], 2, mean, na.rm=TRUE),4)
      tmp[[h]]$eval[[4]][5,] = round(apply(tmp[[h]]$est.MI$simple[,ubcol] - tmp[[h]]$est.MI$simple[,lbcol], 2, mean, na.rm=TRUE),4)
      tmp[[h]]$eval[[5]][5,] = round(apply(tmp[[h]]$est.MI$coinflip[,ubcol] - tmp[[h]]$est.MI$coinflip[,lbcol], 2, mean, na.rm=TRUE),4)
      tmp[[h]]$eval[[6]][5,] = round(apply(tmp[[h]]$est.MI$adaptive[,ubcol] - tmp[[h]]$est.MI$adaptive[,lbcol], 2, mean, na.rm=TRUE),4)
      }
    
    for (h in d3) {     #h: n
      nh = n[h]
      
      # point estimates
      val <- vector(length=6)
      val[1] <- mean(tmp[[h]]$est.com$AUC.hat, na.rm = T)
      val[2] <- mean(tmp[[h]]$est.MI$pmm$AUC.hat, na.rm = T)
      val[3] <- mean(tmp[[h]]$est.MI$logreg$AUC.hat, na.rm = T)
      val[4] <- mean(tmp[[h]]$est.MI$simple$AUC.hat, na.rm = T)
      val[5] <- mean(tmp[[h]]$est.MI$coinflip$AUC.hat, na.rm = T)
      val[6] <- mean(tmp[[h]]$est.MI$adaptive$AUC.hat, na.rm = T)
      
      for (l in d4) {   #l: mi methods
        mh = mi.names[l]
        row.range = index:(index+increment-1)
        col.range = 7:len.col
        result[row.range,1:6] = matrix(c(rep(c(th,ph,rh,nh,mh),each=length(msmt)),msmt),length(msmt))
        result[row.range,col.range] = tmp[[h]]$eval[[l]]
        index=index+increment  # for next
        
        pointest[mmmm,] = c(th,ph,rh,nh,mh, val[l])
        mmmm <- mmmm + 1
        
      }
    }
  }
}
result$n <- as.numeric(result$n)
head(result)
result[result$measure=="CIL" & result$Bm >=1,]

#point estimates for bias check
pointest$MI <- as.factor(pointest$MI)
pointest$MI = factor(pointest$MI,levels(pointest$MI)[c(3,5,4,1,6,2)])
levels(pointest$MI) <- c("complete", "PMM", "LR", "NORM", "NORM(simple)", "NORM(coinflip)")
pointest$AUC.hat <- as.numeric(pointest$AUC.hat)
head(pointest)
pointest.avg <- aggregate(AUC.hat ~ AUC + MI, mean, data=pointest)
ggplot(pointest.avg[pointest.avg$MI %in% c("complete", "PMM", "LR", "NORM"),], aes(MI, AUC.hat)) + geom_point(aes(colour = factor(MI)), size = 4) + facet_grid(. ~ AUC, labeller=label_both) + 
  geom_abline(intercept = 0.8, slope=0) + geom_abline(intercept = 0.9, slope=0) + 
  geom_abline(intercept = 0.95, slope=0) + geom_abline(intercept = 0.99, slope=0) + 
  theme(axis.text.x = element_blank(),axis.ticks = element_blank()) +
  ylab("Average AUC estimate")

## reshaping: wide to long
rst = melt(result, id =c("theta","phi","rho","n","MI","measure"))
names(rst)[7] = "CI.method"
rst$CI.method <- as.factor(rst$CI.method)
rst <- dcast(rst, phi + rho + n + MI + measure + CI.method ~ theta, mean)
names(rst)[7:10] <- paste0("v.",names(rst)[7:10])

## aggregating across phi, rho, n:  table(CI,MI,measure x theta)
rst.avg <- aggregate(cbind(v.0.8,v.0.9,v.0.95,v.0.99) ~ CI.method + MI + measure, mean, data=rst)
aggregate(cbind(v.0.8,v.0.9,v.0.95,v.0.99) ~ CI.method + MI + measure + phi, mean, data=rst)


# Calculate Mean Absolute Error for CP only and combine with rst.avg
rst.MAE <- aggregate(cbind(v.0.8,v.0.9,v.0.95,v.0.99) ~ CI.method + MI + measure, MAE, data=rst[rst$measure=="CP",])
rst.MAE$measure <- "CP.MAE"
rst.avg <- rbind(rst.avg,rst.MAE)

# rst2: Wald, MIset
rst2.avg <- rst.avg[(rst.avg$CI.method %in% Wald) & (rst.avg$MI %in% MIset),]
rst3.avg <- lapply(measurement2, function(x) {a<-rst2.avg[rst2.avg$measure==x,-3]; a})
names(rst3.avg) <- measurement2
rst3.avg <- do.call(cbind,rst3.avg)
rst3.avg <- rst3.avg[,-grep("CI.method",names(rst3.avg))[-1]]
rst3.avg <- rst3.avg[,-grep("MI",names(rst3.avg))[-1]]
# cleansing
rst3.avg[,1:2] <- rst3.avg[,c(2,1)]
names(rst3.avg)[1:2] <- c("MI","CI.method")
rst3.avg[,-(1:2)] <- round(rst3.avg[,-(1:2)],3)

if (FALSE) {   # skip...
  # collapsing three theta's into a column
  for (msmt in measurement2) {      #measurement2 : CP, CP.MAE, LNCP, RNCP, CIL
    if (msmt == "CP") {col=3:6}
    else {col = grep(msmt, names(rst3.avg))}
    rst3.avg[,msmt] <- do.call(paste, c(rst3.avg[col], sep=", "))
  }
  rst3.avg <- rst3.avg[,c("MI","CI.method",measurement2)]
}

for (MI in MIset) {
  print(xtable(rst3.avg[rst3.avg==MI,-1], caption=MI), include.rownames = FALSE)
}

# transpose : wide to long
{
  rst4.avg = t(as.matrix(rst3.avg))
  rst4.avg = cbind (rownames(rst4.avg),rownames(rst4.avg),rst4.avg)
  rst4.avg[,1] = c("MI","CI.method",rep("CP",4),rep("CP.MAE",4),rep("LNCP",4),rep("RNCP",4),rep("CIL",4))
  rst4.avg[,2] = c("MI","CI.method",rep(c(.8,.9,.95,.99),5))
  rownames(rst4.avg) <- NULL
  
  for (MI in MIset) {
    print(xtable(rst4.avg[-1,c(1,2,which(rst4.avg[1,]==MI))], caption=MI), include.rownames = FALSE)
  }  
}



if(FALSE){
  rst <- reshape(rst, v.names = "value", idvar = c("phi","rho","n", "MI", "CI.method"), timevar = "theta", direction = "wide")
 # rst2: Wald, MIset
 rst2 <- rst[(rst$CI.method %in% Wald) & (rst$MI %in% MIset),]
 aggregate(cbind(0.8, 0.9, 0.95, 0.99) ~ CI.method + MI, FUN=MAE, data=rst2[rst2$measure=="CP",])
 
 aggregate(data=rst2[rst2$measure=="CP",], 
           by = data[,c("MI","CI.method","measure")],
           FUN=mean, na.rm=TRUE)
}



## calculating deviations(Actual CP - Nominal CP)
result.dev = melt(result, id =c("theta","phi","rho","n","MI","measure"))
names(result.dev)[7] = "CI.method"
result.dev$CI.method <- as.factor(result.dev$CI.method)
result.dev = result.dev[result.dev$measure=="CP",]
result.dev[,8] <- result.dev[,8]-.95
result.dev[,"measure"] <- "CP.dev"

# Filtering Wald type CI's only
result.Wald <- result.dev[result.dev$CI.method %in% Wald,]




## MAE(mean absolute error)
agg.CMn <- aggregate(value~CI.method+theta+MI+n, FUN=MAE,data=result.Wald)


agg.CM <- aggregate(value~CI.method+theta+MI, FUN=MAE,data=result.Wald)
agg.CMn <- aggregate(value~CI.method+theta+MI+n, FUN=MAE,data=result.Wald)
for (i in 1:6) { for (j in 1:3) {
  with(agg.CMn[agg.CMn$MI==mi.names[i] & agg.CMn$n==n[j] & agg.CMn$CI.method=="Bm",], 
       plot(value ~ theta, type="l", ylab="MAD", ylim=c(0,.2), main=paste0(mi.names[i], ", n= ", n[j])))
  for(l in Wald) {
    with(agg.CMn[agg.CMn$MI==mi.names[i] & agg.CMn$n==n[j] & agg.CMn$CI.method==l,], 
         points(value ~ theta, type="l", col=which(l==Wald)))    
    abline(h=c(.02,.03,.04))
  }
  legend(.85, .2, Wald, pch=16, cex=.5, col=1:6)
}}

plot(result.Wald$value ~ result.Wald$theta, type="l", ylim=c(.9,1),  xlab="theta=.8 .9 .95",)
for(l in mainplayers) {points(b[,CI.methods[l]], type="l", pch=16, col=l)}
abline(h=.95, lty=3)

Wald2 <- Wald[c(4,3,2,5,1)]
a <- result[result$measure %in% c("CP") & result$MI %in% c("pmm", "logreg"),c(names(result)[1:6],Wald2)]
for(i in seq(1,31,by=2))  {plot(1:5, a[i,7:11], ylim=c(.9,1), main=paste(a[i,1:4], collapse="-"));points(1:5, a[i+1,7:11], col="red");abline(h=.95)}
dim(a)




## MSE(mean squared error)
aggregate(value~CI.method, FUN=MSE,data=result.Wald)

## to see kind of bias = mean(nominal - actual CP)
aggregate(value~CI.method, FUN=mean,data=result.Wald)


## rst5
rst5 <- rst[(rst$CI.method %in% Wald) & (rst$MI %in% MIset),]
rst5$MI <- as.factor(rst5$MI)
rst5$MI = factor(rst5$MI,levels(rst5$MI)[c(2,4,3,1)])

# reshaping: wide to long
rst5 <- reshape(rst5, varying=c("v.0.8","v.0.9","v.0.95","v.0.99"), direction="long", idvar=c("phi","rho","n","MI","CI.method", "measure"), sep=".0")
rownames(rst5) <- NULL
names(rst5)[7] = "theta"
names(rst5)[8] = "value"
levels(rst5$MI) <- c("complete", "PMM", "LR", "NORM")
# rst5[,9] <- paste(rst5$MI, rst5$CI.method, sep="-")
# names(rst5)[9] = "method"


rst5[, c("phi", "rho", "theta", "MI", "CI.method", "n")] <- 
  lapply(rst5[, c("phi", "rho", "theta", "MI", "CI.method", "n")], as.factor)
rst5.phrh <- aggregate(value~phi+rho+theta+MI+CI.method+measure, FUN=mean,data=rst5) # across n
rst5.n <- aggregate(value~n+theta+MI+CI.method+measure, FUN=mean,data=rst5) # across phi and rho

p.CP <- ggplot(rst5.phrh[rst5.phrh$measure == "CP",], aes(theta, value, group = CI.method)) + 
  geom_line(aes(color=CI.method)) +
  geom_abline(intercept = .95, slope = 0) +
  geom_point(aes(shape=CI.method, color=CI.method)) + 
  facet_grid(phi+rho ~ MI, labeller=label_both, scales="free_y", space="free_y") +
  # ylim(0.7,1) +
  xlab("AUC") +
  ylab("CP") +
  # ggtitle("The average coverage probability by each method") +
  scale_y_continuous(minor_breaks = seq(0.5, 1, 0.05)) +
  theme_bw()  +
  theme(legend.position="bottom")
  # geom_text(aes(theta, CP, label=paste("phi & rho", labs), group=NULL), size = 4, color = "grey50", data=dat, parse = T)
p.CP


p.CIL <- ggplot(rst5.phrh[rst5.phrh$measure == "CIL",], aes(theta, value, group = CI.method)) + 
  geom_line(aes(color=CI.method)) +
  # geom_abline(intercept = .95, slope = 0) +
  geom_point(aes(shape=CI.method, color=CI.method)) + 
  # facet_grid(phi+rho ~ MI, labeller=label_both, scales="free_y", space="free_y") +
  facet_grid(phi+rho ~ MI, labeller=label_both) +
  xlab("AUC") +
  ylab("CIL") +
  # ggtitle("The average coverage probability by each method") +
  scale_y_continuous(minor_breaks = seq(0.5, 1, 0.05)) +
  theme_bw()  +
  theme(legend.position="bottom")
p.CIL


p.CP.n <- ggplot(rst5.n[rst5.n$measure == "CP",], aes(theta, value, group = CI.method)) + 
  geom_line(aes(color=CI.method)) +
  geom_abline(intercept = .95, slope = 0) +
  geom_point(aes(shape=CI.method, color=CI.method)) + 
  facet_grid(n ~ MI, labeller=label_both, scales="free_y", space="free_y") +
  # ylim(0.7,1) +
  xlab("AUC") +
  ylab("CP") +
  # ggtitle("The average coverage probability by each method") +
  scale_y_continuous(minor_breaks = seq(0.5, 1, 0.05)) +
  theme_bw()  +
  theme(legend.position="bottom")
# geom_text(aes(theta, CP, label=paste("phi & rho", labs), group=NULL), size = 4, color = "grey50", data=dat, parse = T)
p.CP.n


p.CIL.n <- ggplot(rst5.n[rst5.n$measure == "CIL",], aes(theta, value, group = CI.method)) + 
  geom_line(aes(color=CI.method)) +
  # geom_abline(intercept = .95, slope = 0) +
  geom_point(aes(shape=CI.method, color=CI.method)) + 
  # facet_grid(n ~ MI, labeller=label_both, scales="free_y", space="free_y") +
  facet_grid(n ~ MI, labeller=label_both) +
  xlab("AUC") +
  ylab("CIL") +
  # ggtitle("The average coverage probability by each method") +
  scale_y_continuous(minor_breaks = seq(0.5, 1, 0.05)) +
  theme_bw()  +
  theme(legend.position="bottom")
p.CIL.n