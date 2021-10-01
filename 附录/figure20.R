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
BU2 <- read.table("BU2.dat",header = TRUE)
#BU2 <- BU2[BU2$X==1,]
BU2 <- BU2[BU2$X==0,]
n <- nrow(BU2)

FFF <- (BU2[BU2$censor==1,1]+BU2[BU2$censor==1,2])/2
BU2[BU2$censor==1,1] <- FFF
BU2[BU2$censor==1,2] <- FFF
BU2[BU2$censori==1,1] <- 0
BU2[BU2$censori==3,2] <- Inf
FFF <- seq(from=0,to=20,length.out=10000)
U_BU2 <- sort(unique(c(BU2$L,BU2$R[BU2$R<Inf])))
fff <- Turnbull(FFF,BU2[,1],BU2[,2],!BU2[,3]) #cure:0.3357252

non <- Turnbull(U_BU2,BU2[,1],BU2[,2],!BU2[,3])

tb2 <- non$surv
phihat <- 1-tb2[length(U_BU2)]
dat2 <- data.frame(left=BU2[,1],right=BU2[,2],delta1=(BU2[,4]==1),delta2=(BU2[,4]==2),delta3=(BU2[,4]==3))
#===============================================================================
ml_p1 <- (function(theta){f("weibull",dat2,theta,tb2)}) %>%
  maxNM(., start=c(.5,.5), constraints=list(ineqA= matrix(c(1,0,0,1),2,2), ineqB=rep(0,2)))
ml_p2 <- (function(theta){f("gompertz",dat2,theta,tb2)}) %>%
  maxNM(., start=c(0.1,0.1), constraints=list(ineqA= matrix(c(1,0,0,1),2,2), ineqB=rep(0,2)))
ml_p3 <- (function(theta){f("lognormal",dat2,theta,tb2)}) %>%
  maxNM(., start=c(1,2), constraints=list(ineqA= matrix(c(1,0,0,1),2,2), ineqB=rep(0,2)))
ml_p4 <- (function(theta){f("uniform",dat2,theta,tb2)}) %>%
  maxNM(., start=c(0.1,max(dat2$right[dat2$right<Inf])+0.1), constraints=list(ineqA= matrix(c(1,0,0,1),2,2), ineqB=rep(0,2)))

setEPS()
#postscript(file="C:/Users/Administrator/Documents/Statistics/待议论文/cure/GOF/GOF for Turnbull/tex/Fig20.eps")
postscript(file="C:/Users/Administrator/Documents/Statistics/待议论文/cure/GOF/GOF for Turnbull/tex/Fig21.eps")

par(mfrow=c(2,2))
#1,1
plot(FFF,fff$surv,xlim=c(0,20),ylim=c(0,1),type='l',
     xaxs="i", yaxs="i",ylab='Survival',xlab='Time')
y1<-1-phihat+phihat*pweibull(FFF,ml_p1$estimate[1],ml_p1$estimate[2],lower.tail = FALSE)
lines(FFF,y1,xlim=c(0,20),ylim=c(0,1),type='l',
      xaxs="i", yaxs="i",lty=2)
ml_p11 <- (function(theta){f_nc("weibull",dat2,theta)}) %>%
  maxNM(., start=c(.5,.5), constraints=list(ineqA= matrix(c(1,0,0,1),2,2), ineqB=rep(0,2)))
y11 <- pweibull(FFF,ml_p11$estimate[1],ml_p11$estimate[2],lower.tail = FALSE)
lines(FFF,y11,xlim=c(0,20),ylim=c(0,1),type='l',
      xaxs="i", yaxs="i",lty=3)
legend("topright",legend=c("Turnbull","Weibull without cure","Weibull with cure"),lwd=1,lty=1:3)
#1,2
plot(FFF,fff$surv,xlim=c(0,20),ylim=c(0,1),type='l',
     xaxs="i", yaxs="i",ylab='Survival',xlab='Time')
y2<-1-phihat+phihat*exp(-ml_p2$estimate[1]/ml_p2$estimate[2]*(exp(ml_p2$estimate[2]*FFF)-1))
lines(FFF,y2,xlim=c(0,20),ylim=c(0,1),type='l',
      xaxs="i", yaxs="i",lty=2)
ml_p21 <- (function(theta){f_nc("gompertz",dat2,theta)}) %>%
  maxNM(., start=c(0.1,0.1), constraints=list(ineqA= matrix(c(1,0,0,1),2,2), ineqB=rep(0,2)))
y21 <- pgpz(FFF,ml_p21$estimate[1],ml_p21$estimate[2],lower.tail = FALSE)
lines(FFF,y21,xlim=c(0,20),ylim=c(0,1),type='l',
      xaxs="i", yaxs="i",lty=3)
legend("topright",legend=c("Turnbull","Gompertz without cure","Gompertz with cure"),lwd=1,lty=1:3)

#2,1
plot(FFF,fff$surv,xlim=c(0,20),ylim=c(0,1),type='l',
     xaxs="i", yaxs="i",ylab='Survival',xlab='Time')
y3<-1-phihat+phihat*plnorm(FFF,ml_p3$estimate[1],ml_p3$estimate[2],lower.tail = FALSE)
lines(FFF,y3,xlim=c(0,20),ylim=c(0,1),type='l',
      xaxs="i", yaxs="i",lty=2)
ml_p31 <- (function(theta){f_nc("lognormal",dat2,theta)}) %>%
  maxNM(., start=c(1,2), constraints=list(ineqA= matrix(c(1,0,0,1),2,2), ineqB=rep(0,2)))

y31 <- plnorm(FFF,ml_p31$estimate[1],ml_p31$estimate[2],lower.tail = FALSE)
lines(FFF,y31,xlim=c(0,20),ylim=c(0,1),type='l',
      xaxs="i", yaxs="i",lty=3)
legend("topright",legend=c("Turnbull","lognormal without cure","lognormal with cure"),lwd=1,lty=1:3)
#2,2
plot(FFF,fff$surv,xlim=c(0,20),ylim=c(0,1),type='l',
     xaxs="i", yaxs="i",ylab='Survival',xlab='Time')
y4<-1-phihat+phihat*punif(FFF,ml_p4$estimate[1],ml_p4$estimate[2],lower.tail = FALSE)
lines(FFF,y4,xlim=c(0,20),ylim=c(0,1),type='l',
      xaxs="i", yaxs="i",lty=2)
ml_p41 <- (function(theta){f_nc("uniform",dat2,theta)}) %>%
  maxNM(., start=c(0.1,max(dat2$right[dat2$right<Inf])+10), constraints=list(ineqA= matrix(c(1,0,0,1),2,2), ineqB=rep(0,2)))

y41 <- punif(FFF,ml_p41$estimate[1],ml_p41$estimate[2],lower.tail = FALSE)
lines(FFF,y41,xlim=c(0,20),ylim=c(0,1),type='l',
      xaxs="i", yaxs="i",lty=3)
legend("topright",legend=c("Turnbull","uniform without cure","uniform with cure"),lwd=1,lty=1:3)
dev.off()