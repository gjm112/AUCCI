## library
library(reshape)
MAE = function(CP.dev) {mean(abs(CP.dev))}
MSE = function(CP.dev) {sqrt(sum(CP.dev^2))/(length(CP.dev)-1)}
mainplayers = c("Bm", "HM1", "HM2", "NS1", "NW", "NS2", "Mee", "DL", "RG")   #or mainplayers=1:12

## combining evaluations in a sheet
msmt = c("CP","LNCP","RNCP","NaN","CIL", "ZWI")
mi.names=c("complete","pmm","logreg")
n = c(200,100,50)
d1 = 1:length(param1$alphabet$theta); d2 = 1:dim(param1$gamma)[1]; d3 = 1:length(n); d4 = 1:length(mi.names)


# data formatting
len.col = length(CI.methods) + 6   # 5 for th, ph, rh, n, and mi
len.row = length(msmt)*max(d1)*max(d2)*max(d3)*max(d4)
agreg.table = as.data.frame(matrix(NA,len.row,len.col))
names(agreg.table) = c("theta","phi","rho","n","MI","measure",CI.methods)
rownames(agreg.table) = NULL

index=1; increment=length(msmt)
for (i in d1) {         #i: theta*phi
  th = round(param1$alphabet$theta[i],4); ph = param1$alphabet$phi[i]
  for (j in d2) {       #j: rho
    rh = param1$gamma$rho[j]
    tmp <- readRDS(paste0("Simdata_Nov20/sim_data-Nov20-",i,j,".rds"))
    for (h in d3) {     #h: n
      nh = n[h]
      for (l in d4) {   #l: mi methods
        mh = mi.names[l]
        row.range = index:(index+increment-1)
        col.range = 7:len.col
        agreg.table[row.range,1:6] =  matrix(c(rep(c(th,ph,rh,nh,mh),each=length(msmt)),msmt),length(msmt))
        agreg.table[row.range,col.range] = tmp[[h]]$eval[[l]]
        index=index+increment  # for next
      } 
    }
  }
}

## calculating deviances(Actual CP - Nominal CP)
agreg.table.dev = melt(agreg.table, id =c("theta","phi","rho","n","MI","measure"))
names(agreg.table.dev)[7] = "CI.method"
agreg.table.dev = agreg.table.dev[agreg.table.dev$measure=="CP",]
agreg.table.dev[,8] <- agreg.table.dev[,8]-.95
agreg.table.dev[,"measure"] <- "CP.dev"

## MAE(mean absolute error)
aggregate(value~CI.method, FUN=MAE,data=agreg.table.dev)
aggregate(value~CI.method+MI, FUN=MAE,data=agreg.table.dev)
aggregate(value~CI.method+phi, FUN=MAE,data=agreg.table.dev)

## MSE(mean squared error)
aggregate(value~CI.method, FUN=MSE,data=agreg.table.dev)

## to see kind of bias = mean(nominal - actual CP)
aggregate(value~CI.method, FUN=mean,data=agreg.table.dev)


#### plot TBD
plot(agreg.table.dev$value ~ agreg.table.dev$theta, type="l", ylim=c(.9,1), ylab="CP", xlab="theta=.8 .9 .95",main=paste0(mi.names[mi], " phi= ", thph$phi[i], ", rho= ", rho[j], ", n= ", n[k]))
for(l in mainplayers) {points(b[,CI.methods[l]], type="l", pch=16, col=l)}
abline(h=.95, lty=3)
legend(2.5, 1, CI.methods[mainplayers], pch=16, cex=.5, col=mainplayers)  
cat(paste("n=",n[k]),"\n");print(b)


rm(msmt,th, ph, rh, nh, mh, len.col,len.row,index,increment)



ph=1; j=2; mi=1    #ph=0; j:rho, k:n, mi:1comp~3logreg

for (mi in 1:3) {   #mi for complete/pmm/logreg
  cat ("################ dataset: ", mi.names[mi], "################## \n")
  for (k in 1:3) {   #k for n
    cat ("### dataset: ", mi.names[mi],", n=", n[k], "### \n")
    for (i in (c(1:3)+3*ph)) {  # i for theta*phi
      a <- readRDS(paste0("Simdata_Nov19/sim_data-Nov19-",i,j,".rds"))
      b[i -3*(i>3),] <-a[[k]]$eval[[mi]][1,]
    }
    plot(b[,"Bm"], type="l", ylim=c(.9,1), ylab="CP", xlab="theta=.8 .9 .95",main=paste0(mi.names[mi], " phi= ", thph$phi[i], ", rho= ", rho[j], ", n= ", n[k]))
    for(l in mainplayers) {points(b[,CI.methods[l]], type="l", pch=16, col=l)}
    abline(h=.95, lty=3)
    legend(2.5, 1, CI.methods[mainplayers], pch=16, cex=.5, col=mainplayers)  
    cat(paste("n=",n[k]),"\n");print(b)
  }
}