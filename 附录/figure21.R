library(ReIns)
library(maxLik)
library(magrittr)
library(survival)
library(pracma)
library(icensBKL)
source("gpz.R")
source("f.R")
source("f_nocure.R")
#================================================================================
data(breastCancer)
bC <- breastCancer[,-3]
bC <- bC[-42,]
rm(breastCancer)
n <- nrow(bC)

bC$delta1 <- is.na(bC$low)
bC$delta2 <- (!is.na(bC$upp))&(!is.na(bC$low))
bC$delta3 <- is.na(bC$upp)
bC[bC$delta2==1,]$delta2 <- (bC[bC$delta2==1,2]-bC[bC$delta2==1,1])>=5
bC[!rowSums(bC[,3:5]),1] <- bC[!rowSums(bC[,3:5]),2] <- (bC[!rowSums(bC[,3:5]),1]+bC[!rowSums(bC[,3:5]),2])/2
bC[bC$delta1==1,1] <- 0
bC[bC$delta3==1,2] <- Inf

FF <- seq(from=0,to=60,length.out=10000)
U_bc <- sort(unique(c(as.matrix(bC)[,1:2]))[unique(c(as.matrix(bC)[,1:2]))<Inf])
non <- Turnbull(U_bc,bC[,1],bC[,2],rowSums(bC[,3:5]))
tb3 <- non$surv
ff <- Turnbull(FF,bC$low,bC$upp,rowSums(bC[,3:5]))
tb4 <- ff$surv
phihat <- 1-tb3[length(tb3)]
dat3 <- data.frame(left=bC$low,right=bC$upp,delta1=bC$delta1,delta2=bC$delta2,delta3=bC$delta3)
#==========================================================================
ml_p1 <- (function(theta){f("weibull",dat3,theta,tb3)}) %>%
  maxLik(.,start=c(0.1,0.1),method = "NM")
ml_p2 <- (function(theta){f("gompertz",dat3,theta,tb3)}) %>%
  maxLik(.,start=c(0.01,0.01),method = "NM")
ml_p3 <- (function(theta){f("lognormal",dat3,theta,tb3)}) %>%
  maxLik(.,start=c(0,1),method = "NM")
ml_p4 <- (function(theta){f("uniform",dat3,theta,tb3)}) %>%
  maxLik(.,start=c(0,max(dat3$right[dat3$right<Inf])),method = "NM")


setEPS()
postscript(file="C:/Users/Administrator/Documents/Statistics/´ýÒéÂÛÎÄ/cure/GOF/GOF for Turnbull/tex/Fig21.eps")
par(mfrow=c(2,2))


#1,1
plot(FF,tb4,xlim=c(0,60),ylim=c(0,1),type='l',
     xaxs="i", yaxs="i",ylab='Survival',xlab='Time',lty=1)
y1<-1-phihat+phihat*pweibull(FF,ml_p1$estimate[1],ml_p1$estimate[2],lower.tail = FALSE)
lines(FF,y1,xlim=c(0,60),ylim=c(0,1),type='l',
      xaxs="i", yaxs="i",lty=2)
ml_p11 <- (function(theta){f_nc("weibull",dat3,theta)}) %>%
  maxNM(., start=c(.5,.5), constraints=list(ineqA= matrix(c(1,0,0,1),2,2), ineqB=rep(0,2)))
y11 <- pweibull(FF,ml_p11$estimate[1],ml_p11$estimate[2],lower.tail = FALSE)
lines(FF,y11,xlim=c(0,60),ylim=c(0,1),type='l',
      xaxs="i", yaxs="i",lty=3)
legend("topright",legend=c("Turnbull","Weibull without cure","Weibull with cure"),lwd=1,lty=1:3)
#1,2
plot(FF,tb4,xlim=c(0,60),ylim=c(0,1),type='l',
     xaxs="i", yaxs="i",ylab='Survival',xlab='Time',lty=1)
y2<-1-phihat+phihat*exp(-ml_p2$estimate[1]/ml_p2$estimate[2]*(exp(ml_p2$estimate[2]*FF)-1))
lines(FF,y2,xlim=c(0,60),ylim=c(0,1),type='l',
      xaxs="i", yaxs="i",lty=2)
ml_p21 <- (function(theta){f_nc("gompertz",dat3,theta)}) %>%
  maxNM(., start=c(0.1,0.1), constraints=list(ineqA= matrix(c(1,0,0,1),2,2), ineqB=rep(0,2)))
y21 <- pgpz(FF,ml_p21$estimate[1],ml_p21$estimate[2],lower.tail = FALSE)
lines(FF,y21,xlim=c(0,60),ylim=c(0,1),type='l',
      xaxs="i", yaxs="i",lty=3)
legend("topright",legend=c("Turnbull","Gompertz without cure","Gompertz with cure"),lwd=1,lty=1:3)

#2,1
plot(FF,tb4,xlim=c(0,60),ylim=c(0,1),type='l',
     xaxs="i", yaxs="i",ylab='Survival',xlab='Time',lty=1)
y3<-1-phihat+phihat*plnorm(FF,ml_p3$estimate[1],ml_p3$estimate[2],lower.tail = FALSE)
lines(FF,y3,xlim=c(0,60),ylim=c(0,1),type='l',
      xaxs="i", yaxs="i",lty=2)
ml_p31 <- (function(theta){f_nc("lognormal",dat3,theta)}) %>%
  maxNM(., start=c(1,2), constraints=list(ineqA= matrix(c(1,0,0,1),2,2), ineqB=rep(0,2)))

y31 <- plnorm(FF,ml_p31$estimate[1],ml_p31$estimate[2],lower.tail = FALSE)
lines(FF,y31,xlim=c(0,60),ylim=c(0,1),type='l',
      xaxs="i", yaxs="i",lty=3)
legend("topright",legend=c("Turnbull","lognormal without cure","lognormal with cure"),lwd=1,lty=1:3)
#2,2
plot(FF,tb4,xlim=c(0,60),ylim=c(0,1),type='l',
     xaxs="i", yaxs="i",ylab='Survival',xlab='Time',lty=1)
y4<-1-phihat+phihat*punif(FF,ml_p4$estimate[1],ml_p4$estimate[2],lower.tail = FALSE)
lines(FF,y4,xlim=c(0,60),ylim=c(0,1),type='l',
      xaxs="i", yaxs="i",lty=2)
ml_p41 <- (function(theta){f_nc("uniform",dat3,theta)}) %>%
  maxNM(., start=c(0.1,max(dat3$right[dat3$right<Inf])+0.1), constraints=list(ineqA= matrix(c(1,0,0,1),2,2), ineqB=rep(0,2)))

y41 <- punif(FF,ml_p41$estimate[1],ml_p41$estimate[2],lower.tail = FALSE)
lines(FF,y41,xlim=c(0,60),ylim=c(0,1),type='l',
      xaxs="i", yaxs="i",lty=3)
legend("topright",legend=c("Turnbull","uniform without cure","uniform with cure"),lwd=1,lty=1:3)
dev.off()