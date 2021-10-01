#alpha0.3phi0.8 4-7
#alpha0.1phi0.5 8-11
#=================================================================================
source("ini_value.R")
we_500 <- read.csv("C:/Users/Administrator/Documents/Statistics/待议论文/cure/GOF/GOF for Turnbull/tttttry/500/alpha0.3phi0.8(2021)/we.csv",header = F)
go_500 <- read.csv("C:/Users/Administrator/Documents/Statistics/待议论文/cure/GOF/GOF for Turnbull/tttttry/500/alpha0.3phi0.8(2021)/go.csv",header = F)
lo_500 <- read.csv("C:/Users/Administrator/Documents/Statistics/待议论文/cure/GOF/GOF for Turnbull/tttttry/500/alpha0.3phi0.8(2021)/lo.csv",header = F)
un_500 <- na.omit(read.csv("C:/Users/Administrator/Documents/Statistics/待议论文/cure/GOF/GOF for Turnbull/tttttry/500/alpha0.3phi0.8(2021)/un.csv",header = F))

##Weibull
x <- seq(0,8,1e-3)
theta_hat <- apply(we_500[,16:23],2,mean)
setEPS()
postscript(file="C:/Users/Administrator/Documents/Statistics/待议论文/cure/GOF/GOF for Turnbull/tex/Fig4.eps")
par(mfrow=c(1,3))
#Density
plot(x=x,y=dweibull(x,shape_weib,scale_weib),
     type="l",lty=1,ylab="Density",xlab="Time")
lines(x=x,y=dweibull(x,theta_hat[1],theta_hat[2]),type="l",lty=2)
lines(x=x,y=theta_hat[3]*exp(theta_hat[4]*x)*exp(-theta_hat[3]/theta_hat[4]*(exp(theta_hat[4]*x)-1)),type="l",lty=3)
lines(x=x,y=dlnorm(x,theta_hat[5],theta_hat[6]),type="l",lty=4)
lines(x=x,y=dunif(x,theta_hat[7],theta_hat[8]),type="l",lty=5)
legend("topright",legend=c("true","Weibull","Gompertz","lognormal","uniform"),
       lty=1:5,xjust=1)
#Survival
plot(x=x,y=pweibull(x,shape_weib,scale_weib,lower.tail = F),
     type="l",lty=1,ylab="Survival",xlab="Time")
lines(x=x,y=pweibull(x,theta_hat[1],theta_hat[2],lower.tail = F),type="l",lty=2)
lines(x=x,y=exp(-theta_hat[3]/theta_hat[4]*(exp(theta_hat[4]*x)-1)),type="l",lty=3)
lines(x=x,y=plnorm(x,theta_hat[5],theta_hat[6],lower.tail = F),type="l",lty=4)
lines(x=x,y=punif(x,theta_hat[7],theta_hat[8],lower.tail = F),type="l",lty=5)
legend("topright",legend=c("true","Weibull","Gompertz","lognormal","uniform"),
       lty=1:5,xjust=1)
#Hazard
plot(x=x,y=dweibull(x,shape_weib,scale_weib)/pweibull(x,shape_weib,scale_weib,lower.tail = F),
     type="l",lty=1,ylab="Survival",xlab="Time")
lines(x=x,y=dweibull(x,theta_hat[1],theta_hat[2])/pweibull(x,theta_hat[1],theta_hat[2],lower.tail = F),type="l",lty=2)
lines(x=x,y=theta_hat[3]*exp(theta_hat[4]*x),type="l",lty=3)
lines(x=x,y=dlnorm(x,theta_hat[5],theta_hat[6])/plnorm(x,theta_hat[5],theta_hat[6],lower.tail = F),type="l",lty=4)
lines(x=x,y=dunif(x,theta_hat[7],theta_hat[8])/punif(x,theta_hat[7],theta_hat[8],lower.tail = F),type="l",lty=5)
legend("topright",legend=c("true","Weibull","Gompertz","lognormal","uniform"),
       lty=1:5,xjust=1)
dev.off()
#=================================================================================
##Gompertz
setEPS()
postscript(file="C:/Users/Administrator/Documents/Statistics/待议论文/cure/GOF/GOF for Turnbull/tex/Fig5.eps")
theta_hat <- apply(go_500[,16:23],2,mean)
par(mfrow=c(1,3))
#Density
plot(x=x,y=labda*exp(gama*x)*exp(-labda/gama*(exp(gama*x)-1)),
     type="l",lty=1,ylab="Density",xlab="Time")
lines(x=x,y=dweibull(x,theta_hat[1],theta_hat[2]),type="l",lty=2)
lines(x=x,y=theta_hat[3]*exp(theta_hat[4]*x)*exp(-theta_hat[3]/theta_hat[4]*(exp(theta_hat[4]*x)-1)),type="l",lty=3)
lines(x=x,y=dlnorm(x,theta_hat[5],theta_hat[6]),type="l",lty=4)
lines(x=x,y=dunif(x,theta_hat[7],theta_hat[8]),type="l",lty=5)
legend("topright",legend=c("true","Weibull","Gompertz","lognormal","uniform"),
       lty=1:5,xjust=1)
#Survival
plot(x=x,y=exp(-labda/gama*(exp(gama*x)-1)),
     type="l",lty=1,ylab="Survival",xlab="Time")
lines(x=x,y=pweibull(x,theta_hat[1],theta_hat[2],lower.tail = F),type="l",lty=2)
lines(x=x,y=exp(-theta_hat[3]/theta_hat[4]*(exp(theta_hat[4]*x)-1)),type="l",lty=3)
lines(x=x,y=plnorm(x,theta_hat[5],theta_hat[6],lower.tail = F),type="l",lty=4)
lines(x=x,y=punif(x,theta_hat[7],theta_hat[8],lower.tail = F),type="l",lty=5)
legend("topright",legend=c("true","Weibull","Gompertz","lognormal","uniform"),
       lty=1:5,xjust=1)
#Hazard
plot(x=x,y=labda*exp(gama*x),
     type="l",lty=1,ylab="Survival",xlab="Time")
lines(x=x,y=dweibull(x,theta_hat[1],theta_hat[2])/pweibull(x,theta_hat[1],theta_hat[2],lower.tail = F),type="l",lty=2)
lines(x=x,y=theta_hat[3]*exp(theta_hat[4]*x),type="l",lty=3)
lines(x=x,y=dlnorm(x,theta_hat[5],theta_hat[6])/plnorm(x,theta_hat[5],theta_hat[6],lower.tail = F),type="l",lty=4)
lines(x=x,y=dunif(x,theta_hat[7],theta_hat[8])/punif(x,theta_hat[7],theta_hat[8],lower.tail = F),type="l",lty=5)
legend("topright",legend=c("true","Weibull","Gompertz","lognormal","uniform"),
       lty=1:5,xjust=1)
dev.off()
#=================================================================================
##lognormal
setEPS()
postscript(file="C:/Users/Administrator/Documents/Statistics/待议论文/cure/GOF/GOF for Turnbull/tex/Fig6.eps")
theta_hat <- apply(lo_500[,16:23],2,mean)
par(mfrow=c(1,3))
#Density
plot(x=x,y=dlnorm(x,meanlog,sdlog),
     type="l",lty=1,ylab="Density",xlab="Time")
lines(x=x,y=dweibull(x,theta_hat[1],theta_hat[2]),type="l",lty=2)
lines(x=x,y=theta_hat[3]*exp(theta_hat[4]*x)*exp(-theta_hat[3]/theta_hat[4]*(exp(theta_hat[4]*x)-1)),type="l",lty=3)
lines(x=x,y=dlnorm(x,theta_hat[5],theta_hat[6]),type="l",lty=4)
lines(x=x,y=dunif(x,theta_hat[7],theta_hat[8]),type="l",lty=5)
legend("topright",legend=c("true","Weibull","Gompertz","lognormal","uniform"),
       lty=1:5,xjust=1)
#Survival
plot(x=x,y=plnorm(x,meanlog,sdlog,lower.tail = F),
     type="l",lty=1,ylab="Survival",xlab="Time")
lines(x=x,y=pweibull(x,theta_hat[1],theta_hat[2],lower.tail = F),type="l",lty=2)
lines(x=x,y=exp(-theta_hat[3]/theta_hat[4]*(exp(theta_hat[4]*x)-1)),type="l",lty=3)
lines(x=x,y=plnorm(x,theta_hat[5],theta_hat[6],lower.tail = F),type="l",lty=4)
lines(x=x,y=punif(x,theta_hat[7],theta_hat[8],lower.tail = F),type="l",lty=5)
legend("topright",legend=c("true","Weibull","Gompertz","lognormal","uniform"),
       lty=1:5,xjust=1)
#Hazard
plot(x=x,y=dlnorm(x,meanlog,sdlog)/plnorm(x,meanlog,sdlog,lower.tail = F),
     type="l",lty=1,ylab="Survival",xlab="Time")
lines(x=x,y=dweibull(x,theta_hat[1],theta_hat[2])/pweibull(x,theta_hat[1],theta_hat[2],lower.tail = F),type="l",lty=2)
lines(x=x,y=theta_hat[3]*exp(theta_hat[4]*x),type="l",lty=3)
lines(x=x,y=dlnorm(x,theta_hat[5],theta_hat[6])/plnorm(x,theta_hat[5],theta_hat[6],lower.tail = F),type="l",lty=4)
lines(x=x,y=dunif(x,theta_hat[7],theta_hat[8])/punif(x,theta_hat[7],theta_hat[8],lower.tail = F),type="l",lty=5)
legend("topright",legend=c("true","Weibull","Gompertz","lognormal","uniform"),
       lty=1:5,xjust=1)
dev.off()
#=================================================================================
##uniform
setEPS()
postscript(file="C:/Users/Administrator/Documents/Statistics/待议论文/cure/GOF/GOF for Turnbull/tex/Fig7.eps")
theta_hat <- apply(un_500[,16:23],2,mean)

par(mfrow=c(1,3))
#Density
plot(x=x,y=dunif(x,min_unif,max_unif),
     type="l",lty=1,ylab="Density",xlab="Time")
lines(x=x,y=dweibull(x,theta_hat[1],theta_hat[2]),type="l",lty=2)
lines(x=x,y=theta_hat[3]*exp(theta_hat[4]*x)*exp(-theta_hat[3]/theta_hat[4]*(exp(theta_hat[4]*x)-1)),type="l",lty=3)
lines(x=x,y=dlnorm(x,theta_hat[5],theta_hat[6]),type="l",lty=4)
lines(x=x,y=dunif(x,theta_hat[7],theta_hat[8]),type="l",lty=5)
legend("topright",legend=c("true","Weibull","Gompertz","lognormal","uniform"),
       lty=1:5,xjust=1)
#Survival
plot(x=x,y=punif(x,min_unif,max_unif,lower.tail = F),
     type="l",lty=1,ylab="Survival",xlab="Time")
lines(x=x,y=pweibull(x,theta_hat[1],theta_hat[2],lower.tail = F),type="l",lty=2)
lines(x=x,y=exp(-theta_hat[3]/theta_hat[4]*(exp(theta_hat[4]*x)-1)),type="l",lty=3)
lines(x=x,y=plnorm(x,theta_hat[5],theta_hat[6],lower.tail = F),type="l",lty=4)
lines(x=x,y=punif(x,theta_hat[7],theta_hat[8],lower.tail = F),type="l",lty=5)
legend("topright",legend=c("true","Weibull","Gompertz","lognormal","uniform"),
       lty=1:5,xjust=1)
#Hazard
plot(x=x,y=dunif(x,min_unif,max_unif)/punif(x,min_unif,max_unif,lower.tail = F),
     type="l",lty=1,ylab="Survival",xlab="Time")
lines(x=x,y=dweibull(x,theta_hat[1],theta_hat[2])/pweibull(x,theta_hat[1],theta_hat[2],lower.tail = F),type="l",lty=2)
lines(x=x,y=theta_hat[3]*exp(theta_hat[4]*x),type="l",lty=3)
lines(x=x,y=dlnorm(x,theta_hat[5],theta_hat[6])/plnorm(x,theta_hat[5],theta_hat[6],lower.tail = F),type="l",lty=4)
lines(x=x,y=dunif(x,theta_hat[7],theta_hat[8])/punif(x,theta_hat[7],theta_hat[8],lower.tail = F),type="l",lty=5)
legend("topright",legend=c("true","Weibull","Gompertz","lognormal","uniform"),
       lty=1:5,xjust=1)
dev.off()