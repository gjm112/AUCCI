## plot of CP by theta

# for empty shell b
a <- readRDS(paste0("Simdata_Nov15/sim_data-Nov15-11.rds"))
b <- a[[1]]$eval[[1]][1:3,]
rownames(b) <- c("th80","th90","th95")

thph=data.frame(i=1:6, theta=rep(c(.8,.9,.95),2), phi=rep(c(.5,.7),each=3))
mi.names=c("1.complete","2.pmm","3.logreg")

phi=c(.5,.7);rho=c(.5,.7)
# without 5NW,6NS2,10DB,11DG,12CM
mainplayers = c(1:4,7:9)   #or mainplayers=1:12

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

if (FALSE) { # for practice
  head(a[[1]]$est.com[[1]])
  names(a[[1]])
  head(a[[1]]$est.com[[1]])
  a[[3]]$est.MI$pmm[1:10,]
  a[[3]]$eval[[2]]
  a[[1]]$est.MI$logreg[1:10,]
}

for (mi in 1:3) {    #j:rho, k:n, mi:1comp~3logreg
  cat ("################ dataset: ", mi.names[mi], "################## \n")
  for (i in 1:6) {  # i for theta*phi
    for (j in 1:2) {
      cat("##### theta: ",thph$theta[i],"phi: ",thph$phi[i],"rho: ", rho[j], "##### \n")
      for (k in 1:3) {   #k for n
        a <- readRDS(paste0("Simdata_Nov15/sim_data-Nov15-",i,j,".rds"))
        cat("theta: ",thph$theta[i],"phi: ",thph$phi[i],"rho: ", rho[j],"n=",n[k],"\n")
        print(a[[k]]$eval[[mi]])    
      }
    }
  }
}


## Checking Bias
# AUC.hats for complete data
k=3;
mean(a[[k]]$est.com[[1]], na.rm=TRUE); median(a[[k]]$est.com[[1]], na.rm=TRUE); hist(a[[k]]$est.com[[1]], xlim=c(.7,1), breaks=20)
mean(a[[k]]$est.MI$pmm[[1]], na.rm=TRUE); median(a[[k]]$est.MI$pmm[[1]], na.rm=TRUE); hist(a[[k]]$est.MI$pmm[[1]], xlim=c(.7,1), breaks=30)
mean(a[[k]]$est.MI$logreg[[1]], na.rm=TRUE); median(a[[k]]$est.MI$logreg[[1]], na.rm=TRUE); hist(a[[k]]$est.MI$logreg[[1]], xlim=c(.7,1), breaks=20)
