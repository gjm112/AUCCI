datagenerator.normal <- function(n, theta, mu.x=1, sig.x=1, sig.y=2, phi, fixed.n.x = FALSE){
  mu.y = mu.x + sqrt(sig.x^2 + sig.y^2)*qnorm(theta)
  if (fixed.n.x==TRUE) {n.x = floor(n*phi)} else {
    n.x = sum(rbinom(n,1, phi))
  }
  n.y = n - n.x
  x = rnorm(n.x, mu.x, sig.x)
  y = rnorm(n.y, mu.y, sig.y)
  data.frame(disease = c(rep(0,n.x), rep(1,n.y)), marker = c(x,y))
}

CI.methods = c("Bm", "HM1", "HM2", "NS1", "NW", "NS2", "Mee", "RG")
th=0.9
set.seed(10)
sig = seq(1, 3, by=0.1); 
len = length(sig)
vec <- vector(); mat <- data.frame(theta=rep(th,len), sig.y = rep(NA,len), mean=rep(NA,len), sd=rep(NA,len))
temp.CI <- as.data.frame(matrix(NA,1000,(length(CI.methods)*2+1)))
temp.CI2 <- list()
for (i in 1:len){
  for (j in 1:1000) {
    temp <- datagenerator.normal(40, theta=th, sig.y = sig[i], phi=.5,fixed.n.x=TRUE)
    vec[j] <- AUC( data= temp)
    temp.CI[j,] <- CI.i(data=temp, fun=AUCCI, CI.method=CI.methods, type="landscape2")
    print(c(i,j))
  }
  temp.CI2[[i]] <- CI.evaluator(temp.CI, param = data.frame(theta=th), CI.method = CI.methods, na.rm = TRUE, round=4)
  hist(vec)
  mat[i,] <- data.frame(theta=th, sig.y = sig[i], mean=mean(vec), sd=sd(vec))
}
mat
plot(mat$sig.y,mat$sd)

temp.CI3 <- data.frame(matrix(NA,21,length(CI.methods)))
colnames(temp.CI3) <- CI.methods
rownames(temp.CI3) <- seq(1, 3, by=0.1)
for(k in 1:len) {temp.CI3[k,] <- temp.CI2[[k]][1,]}
temp.CI3
#saveRDS(temp.CI3,"HyperscedasticityTest.rds")
#saveRDS(mat,"HyperscedasticityTest-mat.rds")
