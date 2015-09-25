### alpha0 finder for phi ###

phi = 0.7; alpha0 = log(phi/(1-phi))
alpha0 =1.61107
n.sample = 10000
pb <- txtProgressBar(min=0, max = n.sample, char = paste0("="), style=3)

a <- Sys.time(); phi.hat = rep(NA,n.sample)
set.seed(300)
for (i in 1:n.sample){
  setTxtProgressBar(pb,i)
  temp <- datagenerator(n=1000000, alpha0=alpha0, alpha1=alpha1, beta0=beta0, beta1=beta1, beta2=beta2, beta3=beta3, sig=sig, q1=q1, q2=q2, q3=q3, q4=q4, gamma=gamma, mu.V=mu.V, Sigma=Sigma, option="VD") ; Sys.time()-a
  phi.hat[i] = mean(temp$disease)
  if (i %% floor(n.sample/30) == 0) {print(mean(phi.hat, na.rm=TRUE))}
}
mean(phi.hat)-phi;mean(phi.hat) ; sd(phi.hat); sd(phi.hat)/sqrt(n.sample)
Sys.time() - a
## phi.hat.16110 = phi.hat
