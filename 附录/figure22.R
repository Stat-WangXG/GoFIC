setEPS()
postscript(file="C:/Users/Administrator/Documents/Statistics/´ýÒéÂÛÎÄ/cure/GOF/GOF for Turnbull/tex/Fig22.eps")
par(mfrow=c(1,2))
library(ReIns)
library(maxLik)
library(magrittr)
library(survival)
library(pracma)
library(icensBKL)
#1,2
source("f.R")
BU2 <- read.table("BU2.dat",header = TRUE)
BU2 <- BU2[BU2$X==0,]
n <- nrow(BU2)
FFF <- (BU2[BU2$censor==1,1]+BU2[BU2$censor==1,2])/2
BU2[BU2$censor==1,1] <- FFF
BU2[BU2$censor==1,2] <- FFF
BU2[BU2$censori==1,1] <- 0
BU2[BU2$censori==3,2] <- Inf

U_BU2 <- sort(unique(c(BU2$L,BU2$R[BU2$R<Inf])))
non <- Turnbull(U_BU2,BU2[,1],BU2[,2],!BU2[,3])
tb2 <- non$surv
phihat <- 1-tb2[length(U_BU2)]
dat2 <- data.frame(left=BU2[,1],right=BU2[,2],delta1=(BU2[,4]==1),delta2=(BU2[,4]==2),delta3=(BU2[,4]==3))
#==========================================================================
ml_p1 <- (function(theta){f("weibull",dat2,theta,tb2)}) %>%
  maxNM(., start=c(.5,.5), constraints=list(ineqA= matrix(c(1,0,0,1),2,2), ineqB=rep(0,2)))
ml_p2 <- (function(theta){f("gompertz",dat2,theta,tb2)}) %>%
  maxNM(., start=c(0.1,0.1), constraints=list(ineqA= matrix(c(1,0,0,1),2,2), ineqB=rep(0,2)))
ml_p3 <- (function(theta){f("lognormal",dat2,theta,tb2)}) %>%
  maxNM(., start=c(1,2), constraints=list(ineqA= matrix(c(1,0,0,1),2,2), ineqB=rep(0,2)))
ml_p4 <- (function(theta){f("uniform",dat2,theta,tb2)}) %>%
  maxNM(., start=c(0.1,max(dat2$right[dat2$right<Inf])+0.1), constraints=list(ineqA= matrix(c(1,0,0,1),2,2), ineqB=rep(0,2)))

plot(x = 1-phihat+phihat*pweibull(U_BU2,ml_p1$estimate[1],ml_p1$estimate[2],lower.tail = FALSE), y = tb2, 
     type="p", pch=1, xlab = "Parametric fit", ylab = "Turnbull", xlim = c(0.2,1.1), ylim = c(0.2,1.0))
lines(x = 1-phihat+phihat*exp(-ml_p2$estimate[1]/ml_p2$estimate[2]*(exp(ml_p2$estimate[2]*U_BU2)-1)), y = tb2, 
      type="p", pch=2, xlim = c(0.2,1.1))
lines(x = 1-phihat+phihat*plnorm(U_BU2,ml_p3$estimate[1],ml_p3$estimate[2],lower.tail = FALSE), y = tb2,
      type="p", pch=3, xlim = c(0.2,1.1))
lines(x = 1-phihat+phihat*punif(U_BU2,ml_p4$estimate[1],ml_p4$estimate[2],lower.tail = FALSE), y = tb2,
      type="p", pch=4, xlim = c(0.2,1.1))
lines(x = c(0.2,1.2), y=c(0.2,1.2), type="l")
legend(0.2,1,legend=c("Weibull","Gompertz","lognormal","uniform"),pch=1:4,lty=1,cex=0.8)
#1,1
source("f.R")
BU2 <- read.table("BU2.dat",header = TRUE)
BU2 <- BU2[BU2$X==1,]
n <- nrow(BU2)
FFF <- (BU2[BU2$censor==1,1]+BU2[BU2$censor==1,2])/2
BU2[BU2$censor==1,1] <- FFF
BU2[BU2$censor==1,2] <- FFF
BU2[BU2$censori==1,1] <- 0
BU2[BU2$censori==3,2] <- Inf

U_BU2 <- sort(unique(c(BU2$L,BU2$R[BU2$R<Inf])))
non <- Turnbull(U_BU2,BU2[,1],BU2[,2],!BU2[,3])
tb2 <- non$surv
phihat <- 1-tb2[length(U_BU2)]
dat2 <- data.frame(left=BU2[,1],right=BU2[,2],delta1=(BU2[,4]==1),delta2=(BU2[,4]==2),delta3=(BU2[,4]==3))
#==========================================================================
ml_p1 <- (function(theta){f("weibull",dat2,theta,tb2)}) %>%
  maxNM(., start=c(.5,.5), constraints=list(ineqA= matrix(c(1,0,0,1),2,2), ineqB=rep(0,2)))
ml_p2 <- (function(theta){f("gompertz",dat2,theta,tb2)}) %>%
  maxNM(., start=c(0.1,0.1), constraints=list(ineqA= matrix(c(1,0,0,1),2,2), ineqB=rep(0,2)))
ml_p3 <- (function(theta){f("lognormal",dat2,theta,tb2)}) %>%
  maxNM(., start=c(1,2), constraints=list(ineqA= matrix(c(1,0,0,1),2,2), ineqB=rep(0,2)))
ml_p4 <- (function(theta){f("uniform",dat2,theta,tb2)}) %>%
  maxNM(., start=c(0.1,max(dat2$right[dat2$right<Inf])+0.1), constraints=list(ineqA= matrix(c(1,0,0,1),2,2), ineqB=rep(0,2)))

plot(x = 1-phihat+phihat*pweibull(U_BU2,ml_p1$estimate[1],ml_p1$estimate[2],lower.tail = FALSE), y = tb2, 
     type="p", pch=1, xlab = "Parametric fit", ylab = "Turnbull", xlim = c(0.2,1.1), ylim = c(0.2,1.0))
lines(x = 1-phihat+phihat*exp(-ml_p2$estimate[1]/ml_p2$estimate[2]*(exp(ml_p2$estimate[2]*U_BU2)-1)), y = tb2, 
      type="p", pch=2, xlim = c(0.2,1.1))
lines(x = 1-phihat+phihat*plnorm(U_BU2,ml_p3$estimate[1],ml_p3$estimate[2],lower.tail = FALSE), y = tb2,
      type="p", pch=3, xlim = c(0.2,1.1))
lines(x = 1-phihat+phihat*punif(U_BU2,ml_p4$estimate[1],ml_p4$estimate[2],lower.tail = FALSE), y = tb2,
      type="p", pch=4, xlim = c(0.2,1.1))
lines(x = c(0.2,1.2), y=c(0.2,1.2), type="l")
legend(0.2,1,legend=c("Weibull","Gompertz","lognormal","uniform"),pch=1:4,lty=1,cex=0.8)
#1,2
# data(breastCancer)
# bC <- breastCancer[,-3]
# bC <- bC[-42,]
# rm(breastCancer)
# n <- nrow(bC)
# 
# bC$delta1 <- is.na(bC$low)
# bC$delta2 <- (!is.na(bC$upp))&(!is.na(bC$low))
# bC$delta3 <- is.na(bC$upp)
# bC[bC$delta2==1,]$delta2 <- (bC[bC$delta2==1,2]-bC[bC$delta2==1,1])>=5
# bC[!rowSums(bC[,3:5]),1] <- bC[!rowSums(bC[,3:5]),2] <- (bC[!rowSums(bC[,3:5]),1]+bC[!rowSums(bC[,3:5]),2])/2
# bC[bC$delta1==1,1] <- 0
# bC[bC$delta3==1,2] <- Inf
# 
# U_aid <- sort(unique(c(as.matrix(bC)[,1:2]))[unique(c(as.matrix(bC)[,1:2]))<Inf])
# non <- Turnbull(U_aid,bC[,1],bC[,2],rowSums(bC[,3:5]))
# tb3 <- non$surv
# phihat3 <- 1-tb3[length(tb3)]
# dat3 <- data.frame(left=bC$low,right=bC$upp,delta1=bC$delta1,delta2=bC$delta2,delta3=bC$delta3)
# #==========================================================================
# ml_p1 <- (function(theta){f("weibull",dat3,theta,tb3)}) %>%
#   maxLik(.,start=c(0.1,0.1),method = "NM")
# ml_p2 <- (function(theta){f("gompertz",dat3,theta,tb3)}) %>%
#   maxLik(.,start=c(0.01,0.01),method = "NM")
# ml_p3 <- (function(theta){f("lognormal",dat3,theta,tb3)}) %>%
#   maxLik(.,start=c(0,1),method = "NM")
# ml_p4 <- (function(theta){f("uniform",dat3,theta,tb3)}) %>%
#   maxLik(.,start=c(0,max(dat3$right[dat3$right<Inf])),method = "NM")
# 
# plot(x = 1-phihat3+phihat3*pweibull(U_aid,ml_p1$estimate[1],ml_p1$estimate[2],lower.tail = FALSE), y = tb3, 
#      type="p", pch=1, xlab = "Parametric fit", ylab = "Turnbull", xlim = c(0.2,1.1))
# lines(x = 1-phihat3+phihat3*exp(-ml_p2$estimate[1]/ml_p2$estimate[2]*(exp(ml_p2$estimate[2]*U_aid)-1)), y = tb3, 
#       type="p", pch=2, xlim = c(0.2,1.1))
# lines(x = 1-phihat3+phihat3*plnorm(U_aid,ml_p3$estimate[1],ml_p3$estimate[2],lower.tail = FALSE), y = tb3,
#       type="p", pch=3, xlim = c(0.2,1.1))
# lines(x = 1-phihat3+phihat3*punif(U_aid,ml_p4$estimate[1],ml_p4$estimate[2],lower.tail = FALSE), y = tb3,
#       type="p", pch=4, xlim = c(0.2,1.1))
# lines(x = c(0.2,1.2), y=c(0.2,1.2), type="l")
# legend(0.2,1,legend=c("Weibull","Gompertz","lognormal","uniform"),pch=1:4,lty=1,cex=0.8)
dev.off()