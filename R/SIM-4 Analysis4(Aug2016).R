## library
library(reshape); library(reshape2)
library(xtable);
library(ggplot2)
MAE = function(CP.dev , theta=.95) {mean(abs(CP.dev - theta))}
MSE = function(CP.dev) {sqrt(sum((CP.dev-theta)^2))/(length(CP.dev)-1)}
Wald <- CI.methods
MIset <- c("complete","naive","pmm","logreg", "adaptive")
measurement <- c("CP", "LNCP", "RNCP", "CIL", "NaN", "ZWI")
measurement2 <- c(measurement,"CP.MAE")[c(1,7,2,3,4)]
CI.methods = c("Bm", "HM1", "HM2", "NW", "DL")

## combining evaluations in a sheet
msmt = c("CP","LNCP","RNCP","NaN","CIL", "ZWI")
mi.names=c("complete","naive","pmm","logreg","simple","coinflip","adaptive")
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
    tmp <- readRDS(paste0("R/Simdata5/sim_data-Aug01-",i,j,".rds"))
    # CI lengths to be restricted within [0,1]
    for (h in d3){
      # truncation
      tmp[[h]]$est.com[tmp[[h]]$est.com>1] = 1 ; tmp[[h]]$est.com[tmp[[h]]$est.com<0] = 0
      tmp[[h]]$est.na[tmp[[h]]$est.na>1] = 1 ; tmp[[h]]$est.na[tmp[[h]]$est.na<0] = 0
      tmp[[h]]$est.MI$pmm[tmp[[h]]$est.MI$pmm>1] = 1 ; tmp[[h]]$est.MI$pmm[tmp[[h]]$est.MI$pmm<0] = 0
      tmp[[h]]$est.MI$logreg[tmp[[h]]$est.MI$logreg>1] = 1 ; tmp[[h]]$est.MI$logreg[tmp[[h]]$est.MI$logreg<0] = 0
      tmp[[h]]$est.MI$simple[tmp[[h]]$est.MI$simple>1] = 1 ; tmp[[h]]$est.MI$simple[tmp[[h]]$est.MI$simple<0] = 0
      tmp[[h]]$est.MI$coinflip[tmp[[h]]$est.MI$coinflip>1] = 1 ; tmp[[h]]$est.MI$coinflip[tmp[[h]]$est.MI$coinflip<0] = 0
      tmp[[h]]$est.MI$adaptive[tmp[[h]]$est.MI$adaptive>1] = 1 ; tmp[[h]]$est.MI$adaptive[tmp[[h]]$est.MI$adaptive<0] = 0
      # getting average and then rounding to 4 decimal points
      lgth = length(CI.methods); lbcol = seq(2,2*lgth+1,2); ubcol = seq(3,2*lgth+1,2)
      tmp[[h]]$eval[[1]][5,] = round(apply(tmp[[h]]$est.com[,ubcol] - tmp[[h]]$est.com[,lbcol], 2, mean, na.rm=TRUE),4)
      tmp[[h]]$eval[[2]][5,] = round(apply(tmp[[h]]$est.na[,ubcol] - tmp[[h]]$est.na[,lbcol], 2, mean, na.rm=TRUE),4)
      tmp[[h]]$eval[[3]][5,] = round(apply(tmp[[h]]$est.MI$pmm[,ubcol] - tmp[[h]]$est.MI$pmm[,lbcol], 2, mean, na.rm=TRUE),4)
      tmp[[h]]$eval[[4]][5,] = round(apply(tmp[[h]]$est.MI$logreg[,ubcol] - tmp[[h]]$est.MI$logreg[,lbcol], 2, mean, na.rm=TRUE),4)
      tmp[[h]]$eval[[5]][5,] = round(apply(tmp[[h]]$est.MI$simple[,ubcol] - tmp[[h]]$est.MI$simple[,lbcol], 2, mean, na.rm=TRUE),4)
      tmp[[h]]$eval[[6]][5,] = round(apply(tmp[[h]]$est.MI$coinflip[,ubcol] - tmp[[h]]$est.MI$coinflip[,lbcol], 2, mean, na.rm=TRUE),4)
      tmp[[h]]$eval[[7]][5,] = round(apply(tmp[[h]]$est.MI$adaptive[,ubcol] - tmp[[h]]$est.MI$adaptive[,lbcol], 2, mean, na.rm=TRUE),4)
      }
    
    for (h in d3) {     #h: n
      nh = n[h]
      
      # point estimates
      val <- vector(length=7)
      val[1] <- mean(tmp[[h]]$est.com$AUC.hat, na.rm = T)
      val[2] <- mean(tmp[[h]]$est.na$AUC.hat, na.rm = T)
      val[3] <- mean(tmp[[h]]$est.MI$pmm$AUC.hat, na.rm = T)
      val[4] <- mean(tmp[[h]]$est.MI$logreg$AUC.hat, na.rm = T)
      val[5] <- mean(tmp[[h]]$est.MI$simple$AUC.hat, na.rm = T)
      val[6] <- mean(tmp[[h]]$est.MI$coinflip$AUC.hat, na.rm = T)
      val[7] <- mean(tmp[[h]]$est.MI$adaptive$AUC.hat, na.rm = T)
      
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

#point estimates for bias check
pointest$MI <- as.factor(pointest$MI)
pointest$MI = factor(pointest$MI,levels(pointest$MI)[c(3,5,6,4,1,7,2)])
levels(pointest$MI) <- c("complete","naive", "PMM", "LR", "NORM", "NORM(simple)", "NORM(coinflip)")
pointest$AUC.hat <- as.numeric(pointest$AUC.hat)
head(pointest)
pointest.avg <- aggregate(AUC.hat ~ AUC + MI, mean, data=pointest)

AUC.data <- data.frame(AUC.hat=c(.8,.9,.95,.99), AUC=c(.8,.9,.95,.99))
ggplot(pointest.avg[pointest.avg$MI %in% c("naive","complete", "PMM", "LR", "NORM"),], 
       aes(MI, AUC.hat, label = round(AUC.hat,2))) + 
  geom_text(vjust = 0, nudge_y = 0.005, check_overlap = FALSE)+
  facet_grid(. ~ AUC, labeller=label_both) +
  geom_point(aes(colour = factor(MI), shape=factor(MI)), size = 3) + 
  geom_segment(aes(x=0, xend=5, y=AUC.hat, yend=AUC.hat, label=NULL), size=1,  
               data=AUC.data) +
  theme(axis.text.x = element_blank(),axis.ticks = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "grey80")) +
  scale_x_discrete(NULL) +
  ylab("Average AUC estimate")
ggsave("R/plot_bias.png", width = 200, height = 100, units = "mm")


## reshaping: wide to long
rst = melt(result, id =c("theta","phi","rho","n","MI","measure"))
names(rst)[7] = "CI.method"
rst$CI.method <- as.factor(rst$CI.method)
rst <- dcast(rst, phi + rho + n + MI + measure + CI.method ~ theta, mean)
names(rst)[7:10] <- paste0("v.",names(rst)[7:10])

## aggregating across phi, rho, n:  table(CI,MI,measure x theta)
rst.avg <- aggregate(cbind(v.0.8,v.0.9,v.0.95,v.0.99) ~ CI.method + MI + measure, mean, data=rst)
# aggregate(cbind(v.0.8,v.0.9,v.0.95,v.0.99) ~ CI.method + MI + measure + phi, mean, data=rst)


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

# transpose : wide to long
rst4.avg = t(as.matrix(rst3.avg))
rst4.avg = cbind (rownames(rst4.avg),rownames(rst4.avg),rst4.avg)
rst4.avg[,1] = c("MI","CI.method",rep("CP",4),rep("CP.MAE",4),rep("LNCP",4),rep("RNCP",4),rep("CIL",4))
rst4.avg[,2] = c("MI","CI.method",rep(c(.8,.9,.95,.99),5))
rownames(rst4.avg) <- NULL

for (MI in MIset) {
  print(xtable(rst4.avg[-1,c(1,2,which(rst4.avg[1,]==MI))], caption=MI), include.rownames = FALSE)
}


## rst5
rst5 <- rst[(rst$CI.method %in% Wald) & (rst$MI %in% MIset),]
rst5$MI <- as.factor(rst5$MI)
rst5$MI = factor(rst5$MI,levels(rst5$MI)[c(2,4,5,3,1)])

# reshaping: wide to long
rst5 <- reshape(rst5, varying=c("v.0.8","v.0.9","v.0.95","v.0.99"), direction="long", idvar=c("phi","rho","n","MI","CI.method", "measure"), sep=".0")
rownames(rst5) <- NULL
names(rst5)[7] = "theta"
names(rst5)[8] = "value"
levels(rst5$MI) <- c("complete data", "naive analysis", "PMM", "LR", "NORM")
# rst5[,9] <- paste(rst5$MI, rst5$CI.method, sep="-")
# names(rst5)[9] = "method"


rst5[, c("phi", "rho", "theta", "MI", "CI.method", "n")] <- 
  lapply(rst5[, c("phi", "rho", "theta", "MI", "CI.method", "n")], as.factor)

rst5$phrh[rst5$phi==0.5 & rst5$rho==0.5] = "phi == '.5'~rho == '.5'"
rst5$phrh[rst5$phi==0.5 & rst5$rho==0.7] = "phi == '.5'~rho == '.7'"
rst5$phrh[rst5$phi==0.7 & rst5$rho==0.5] = "phi == '.7'~rho == '.5'"
rst5$phrh[rst5$phi==0.7 & rst5$rho==0.7] = "phi == '.7'~rho == '.7'"

rst5.phrh <- aggregate(value~phrh+theta+MI+CI.method+measure, FUN=mean,data=rst5) # across n
rst5.n <- aggregate(value~n+theta+MI+CI.method+measure, FUN=mean,data=rst5) # across phi and rho
#rst5.phrh$phi2 <- factor(rst5.phrh$phi2, labels = c("phi = 0.5", "phi = 0.7"))
#rst5.phrh$rho2 <- factor(rst5.phrh$rho2, labels = c("rho approx 0.5", "rho approx 0.7"))

p.CP <- ggplot(rst5.phrh[rst5.phrh$measure == "CP",], aes(theta, value, group = CI.method)) + 
  geom_line(aes(color=CI.method, linetype=CI.method)) +
  geom_abline(intercept = .95, slope = 0) +
  geom_point(aes(shape=CI.method, color=CI.method)) + 
  #facet_grid(phi+rho ~ MI, labeller=label_both, scales="free_y", space="free_y") +
  facet_grid(phrh ~ MI, labeller=labeller(.rows = label_parsed, .cols = label_value)) +
  # ylim(0.7,1) +
  xlab("AUC") +
  ylab("CP") +
  # ggtitle("The average coverage probability by each method") +
  scale_y_continuous(minor_breaks = seq(0.5, 1, 0.05)) +
  scale_color_discrete(name="CI methods") +
  scale_shape_discrete(name="CI methods") +
  scale_linetype_discrete(name="CI methods") +
  theme_bw()  +
  theme(legend.position="bottom")
  # geom_text(aes(theta, CP, label=paste("phi & rho", labs), group=NULL), size = 4, color = "grey50", data=dat, parse = T)
p.CP
ggsave("R/plot_CP_phirho.png", width = 200, height = 120, units = "mm")


p.CIL <- ggplot(rst5.phrh[rst5.phrh$measure == "CIL",], aes(theta, value, group = CI.method)) + 
  geom_line(aes(color=CI.method, linetype=CI.method)) +
  # geom_abline(intercept = .95, slope = 0) +
  geom_point(aes(shape=CI.method, color=CI.method)) + 
  # facet_grid(phi+rho ~ MI, labeller=label_both, scales="free_y", space="free_y") +
  facet_grid(phrh ~ MI, labeller=labeller(.rows = label_parsed, .cols = label_value)) +
  xlab("AUC") +
  ylab("CIL") +
  # ggtitle("The average coverage probability by each method") +
  scale_y_continuous(minor_breaks = seq(0.5, 1, 0.05)) +
  scale_color_discrete(name="CI methods") +
  scale_shape_discrete(name="CI methods") +
  scale_linetype_discrete(name="CI methods") +
  theme_bw()  +
  theme(legend.position="bottom")
p.CIL
ggsave("R/plot_CIL_phirho.png", width = 200, height = 120, units = "mm")



p.CP.n <- ggplot(rst5.n[rst5.n$measure == "CP",], aes(theta, value, group = CI.method)) + 
  geom_line(aes(color=CI.method, linetype=CI.method)) +
  geom_abline(intercept = .95, slope = 0) +
  geom_point(aes(shape=CI.method, color=CI.method)) + 
  # facet_grid(n ~ MI, labeller=label_both, scales="free_y", space="free_y") +
  facet_grid(n ~ MI, labeller=labeller(.rows = label_both, .cols = label_value)) +
  # ylim(0.7,1) +
  xlab("AUC") +
  ylab("CP") +
  # ggtitle("The average coverage probability by each method") +
  scale_y_continuous(minor_breaks = seq(0.5, 1, 0.05)) +
  scale_color_discrete(name="CI methods") +
  scale_shape_discrete(name="CI methods") +
  scale_linetype_discrete(name="CI methods") +
  theme_bw()  +
  theme(legend.position="bottom")
# geom_text(aes(theta, CP, label=paste("phi & rho", labs), group=NULL), size = 4, color = "grey50", data=dat, parse = T)
p.CP.n
ggsave("R/plot_CP_n.png", width = 200, height = 100, units = "mm")



p.CIL.n <- ggplot(rst5.n[rst5.n$measure == "CIL",], aes(theta, value, group = CI.method)) + 
  geom_line(aes(color=CI.method, linetype=CI.method)) +
  # geom_abline(intercept = .95, slope = 0) +
  geom_point(aes(shape=CI.method, color=CI.method)) + 
  # facet_grid(n ~ MI, labeller=label_both, scales="free_y", space="free_y") +
  facet_grid(n ~ MI, labeller=labeller(.rows = label_both, .cols = label_value)) +
  xlab("AUC") +
  ylab("CIL") +
  # ggtitle("The average coverage probability by each method") +
  scale_y_continuous(minor_breaks = seq(0.5, 1, 0.05)) +
  scale_color_discrete(name="CI methods") +
  scale_shape_discrete(name="CI methods") +
  scale_linetype_discrete(name="CI methods") +
  theme_bw()  +
  theme(legend.position="bottom")
p.CIL.n
ggsave("R/plot_CIL_n.png", width = 200, height = 100, units = "mm")


p.CP.n.PMMLR <- ggplot(rst5.n[rst5.n$measure == "CP"&(rst5.n$MI == "PMM"|rst5.n$MI == "LR"),], aes(theta, value, group = MI)) + 
  geom_line(aes(color=MI, linetype=MI)) +
  geom_abline(intercept = .95, slope = 0) +
  geom_point(aes(shape=MI, color=MI)) + 
  # facet_grid(n ~ MI, labeller=label_both, scales="free_y", space="free_y") +
  facet_grid(n ~ CI.method, labeller=labeller(.rows = label_both, .cols = label_value)) +
  # ylim(0.7,1) +
  xlab("AUC") +
  ylab("CP") +
  # ggtitle("The average coverage probability by each method") +
  scale_y_continuous(minor_breaks = seq(0.5, 1, 0.05)) +
  scale_color_discrete(name="MI techniques") +
  scale_shape_discrete(name="MI techniques") +
  scale_linetype_discrete(name="MI techniques") +
  theme_bw()  +
  theme(legend.position="bottom")
# geom_text(aes(theta, CP, label=paste("phi & rho", labs), group=NULL), size = 4, color = "grey50", data=dat, parse = T)
p.CP.n.PMMLR
ggsave("R/plot_CP_n_PMMLR.png", width = 200, height = 100, units = "mm")


p.CIL.n.PMMLR <- ggplot(rst5.n[rst5.n$measure == "CIL"&(rst5.n$MI == "PMM"|rst5.n$MI == "LR"),], aes(theta, value, group = MI)) + 
  geom_line(aes(color=MI, linetype=MI)) +
  # geom_abline(intercept = .95, slope = 0) +
  geom_point(aes(shape=MI, color=MI)) + 
  # facet_grid(n ~ MI, labeller=label_both, scales="free_y", space="free_y") +
  facet_grid(n ~ CI.method, labeller=labeller(.rows = label_both, .cols = label_value)) +
  xlab("AUC") +
  ylab("CIL") +
  # ggtitle("The average coverage probability by each method") +
  scale_y_continuous(minor_breaks = seq(0.5, 1, 0.05)) +
  scale_color_discrete(name="MI techniques") +
  scale_shape_discrete(name="MI techniques") +
  scale_linetype_discrete(name="MI techniques") +
  theme_bw()  +
  theme(legend.position="bottom")
p.CIL.n.PMMLR
ggsave("R/plot_CIL_n_PMMLR.png", width = 200, height = 100, units = "mm")

