## library
library(reshape)
MAE = function(CP.dev) {mean(abs(CP.dev))}
MSE = function(CP.dev) {sqrt(sum(CP.dev^2))/(length(CP.dev)-1)}
mainplayers = c("Bm", "HM1", "HM2", "NS1", "NW", "NS2", "Mee", "DL", "RG")   #or mainplayers=1:12

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
for (i in d1) {         #i: theta*phi
  th = round(param1$alphabet$theta[i],4); ph = param1$alphabet$phi[i]
  for (j in d2) {       #j: rho
    rh = param1$gamma$rho[j]
    tmp <- readRDS(paste0("Simdata_Nov30/sim_data-Nov30-",i,j,".rds"))
    for (h in d3) {     #h: n
      nh = n[h]
      for (l in d4) {   #l: mi methods
        mh = mi.names[l]
        row.range = index:(index+increment-1)
        col.range = 7:len.col
        result[row.range,1:6] =  matrix(c(rep(c(th,ph,rh,nh,mh),each=length(msmt)),msmt),length(msmt))
        result[row.range,col.range] = tmp[[h]]$eval[[l]]
        index=index+increment  # for next
      } 
    }
  }
}
result$n <- as.numeric(result$n)

## calculating deviations(Actual CP - Nominal CP)
result.dev = melt(result, id =c("theta","phi","rho","n","MI","measure"))
names(result.dev)[7] = "CI.method"
result.dev$CI.method <- as.factor(result.dev$CI.method)
result.dev = result.dev[result.dev$measure=="CP",]
result.dev[,8] <- result.dev[,8]-.95
result.dev[,"measure"] <- "CP.dev"

# Filtering Wald type CI's only
Wald <- levels(agreg.table.dev$CI.method)[c(1:3,5,8)]
result.Wald <- result.dev[result.dev$CI.method %in% Wald,]

## MAE(mean absolute error)
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

plot(result.Wald$value ~ result.Wald$theta, type="l", ylim=c(.9,1), , xlab="theta=.8 .9 .95",)
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