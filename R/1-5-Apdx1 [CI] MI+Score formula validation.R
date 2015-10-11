dat.gen.p <- function(n,phi, q2, gamma){
  V = rnorm(n)
  D = rbinom(n, 1, phi)
      Q2 = quantile(V, q2)
      rho = (V <= Q2)*gamma
      R = rbinom(n, 1, rho)
      DR = ifelse(R==0, D, NA)    
      return(data.frame(disease=D,V=V, R=R, diseaseR=DR))
}

n1=100; m1=10
temp2= dat.gen.p(n=n1,phi=.9,q2=.8,gamma=.5)
temp2.imp = mice(temp2[,-1],m=m1,printFlag=F)
p.hat.i = rep(NA,m1)              #p.hat(i)
for (i in 1:m1) (p.hat.i[i] = mean(complete(temp2.imp,i)$diseaseR))
V.hat.i = p.hat.i*(1-p.hat.i)/(n1-1)
W.hat = mean(V.hat.i)
B.hat = var(p.hat.i)
T.hat = W.hat + (m1+1)/m1*B.hat
mean(temp2$disease)              #p.hat(comp)
mean(p.hat.i)
W.hat; B.hat; T.hat

# Depending on phi, W varies, and B also varies!!
# for phi near 0.5, W and B are both high
# for phi near 0 or 1, W and B are both very small
# -> so B is dependent on phi (and the degree of missingness of data as well)
