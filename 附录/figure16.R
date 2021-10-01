library(doParallel)
library(snowfall)
#alpha0.1phi0.8 4-7
#alpha0.3phi0.8 8-11
#alpha0.1phi0.5 12-15
#alpha0.3phi0.5 16-19
#=================================================================================
sfInit(parallel = TRUE, cpus = detectCores()-1)
sfSource("allfun.R")
#=================================================================================
##Weibull
x <- seq(0,8,1e-3)
theta_hat <- apply(matrix(unlist(sfLapply(1:tt, bupa)),8,tt),1,mean)
setEPS()
postscript(file="Fig16.eps")
par(mfrow=c(1,3))
#Density
plot(x=x,y=dweibull(x,shape_weib,scale_weib),
     type="l",lty=1,col=1,ylab="Density",xlab="Time")
lines(x=x,y=dweibull(x,theta_hat[1],theta_hat[2]),type="l",lty=1,col=2)
lines(x=x,y=theta_hat[3]*exp(theta_hat[4]*x)*exp(-theta_hat[3]/theta_hat[4]*(exp(theta_hat[4]*x)-1)),type="l",lty=1,col=3)
lines(x=x,y=dlnorm(x,theta_hat[5],theta_hat[6]),type="l",lty=1,col=4)
lines(x=x,y=dunif(x,theta_hat[7],theta_hat[8]),type="l",lty=1,col=7)
legend("topright",legend=c("true","Weibull","Gompertz","lognormal","uniform"),
       lty=1,col=c(1:4,7),xjust=1)
#Survival
plot(x=x,y=pweibull(x,shape_weib,scale_weib,lower.tail = F),
     type="l",lty=1,col=1,ylab="Survival",xlab="Time")
lines(x=x,y=pweibull(x,theta_hat[1],theta_hat[2],lower.tail = F),type="l",lty=1,col=2)
lines(x=x,y=exp(-theta_hat[3]/theta_hat[4]*(exp(theta_hat[4]*x)-1)),type="l",lty=1,col=3)
lines(x=x,y=plnorm(x,theta_hat[5],theta_hat[6],lower.tail = F),type="l",lty=1,col=4)
lines(x=x,y=punif(x,theta_hat[7],theta_hat[8],lower.tail = F),type="l",lty=1,col=7)
legend("topright",legend=c("true","Weibull","Gompertz","lognormal","uniform"),
       lty=1,col=c(1:4,7),xjust=1)
#Hazard
plot(x=x,y=dweibull(x,shape_weib,scale_weib)/pweibull(x,shape_weib,scale_weib,lower.tail = F),
     type="l",lty=1,col=1,ylab="Survival",xlab="Time")
lines(x=x,y=dweibull(x,theta_hat[1],theta_hat[2])/pweibull(x,theta_hat[1],theta_hat[2],lower.tail = F),type="l",lty=1,col=2)
lines(x=x,y=theta_hat[3]*exp(theta_hat[4]*x),type="l",lty=1,col=3)
lines(x=x,y=dlnorm(x,theta_hat[5],theta_hat[6])/plnorm(x,theta_hat[5],theta_hat[6],lower.tail = F),type="l",lty=1,col=4)
lines(x=x,y=dunif(x,theta_hat[7],theta_hat[8])/punif(x,theta_hat[7],theta_hat[8],lower.tail = F),type="l",lty=1,col=7)
legend("topright",legend=c("true","Weibull","Gompertz","lognormal","uniform"),
       lty=1,col=c(1:4,7),xjust=1)
dev.off()
#=================================================================================
##Gompertz
setEPS()
postscript(file="Fig17.eps")
theta_hat <- apply(matrix(unlist(sfLapply(1:tt, bupa2)),8,tt),1,mean)
par(mfrow=c(1,3))
#Density
plot(x=x,y=labda*exp(gama*x)*exp(-labda/gama*(exp(gama*x)-1)),
     type="l",lty=1,col=1,ylab="Density",xlab="Time")
lines(x=x,y=dweibull(x,theta_hat[1],theta_hat[2]),type="l",lty=1,col=2)
lines(x=x,y=theta_hat[3]*exp(theta_hat[4]*x)*exp(-theta_hat[3]/theta_hat[4]*(exp(theta_hat[4]*x)-1)),type="l",lty=1,col=3)
lines(x=x,y=dlnorm(x,theta_hat[5],theta_hat[6]),type="l",lty=1,col=4)
lines(x=x,y=dunif(x,theta_hat[7],theta_hat[8]),type="l",lty=1,col=7)
legend("topright",legend=c("true","Weibull","Gompertz","lognormal","uniform"),
       lty=1,col=c(1:4,7),xjust=1)
#Survival
plot(x=x,y=exp(-labda/gama*(exp(gama*x)-1)),
     type="l",lty=1,col=1,ylab="Survival",xlab="Time")
lines(x=x,y=pweibull(x,theta_hat[1],theta_hat[2],lower.tail = F),type="l",lty=1,col=2)
lines(x=x,y=exp(-theta_hat[3]/theta_hat[4]*(exp(theta_hat[4]*x)-1)),type="l",lty=1,col=3)
lines(x=x,y=plnorm(x,theta_hat[5],theta_hat[6],lower.tail = F),type="l",lty=1,col=4)
lines(x=x,y=punif(x,theta_hat[7],theta_hat[8],lower.tail = F),type="l",lty=1,col=7)
legend("topright",legend=c("true","Weibull","Gompertz","lognormal","uniform"),
       lty=1,col=c(1:4,7),xjust=1)
#Hazard
plot(x=x,y=labda*exp(gama*x),
     type="l",lty=1,col=1,ylab="Survival",xlab="Time")
lines(x=x,y=dweibull(x,theta_hat[1],theta_hat[2])/pweibull(x,theta_hat[1],theta_hat[2],lower.tail = F),type="l",lty=1,col=2)
lines(x=x,y=theta_hat[3]*exp(theta_hat[4]*x),type="l",lty=1,col=3)
lines(x=x,y=dlnorm(x,theta_hat[5],theta_hat[6])/plnorm(x,theta_hat[5],theta_hat[6],lower.tail = F),type="l",lty=1,col=4)
lines(x=x,y=dunif(x,theta_hat[7],theta_hat[8])/punif(x,theta_hat[7],theta_hat[8],lower.tail = F),type="l",lty=1,col=7)
legend("topright",legend=c("true","Weibull","Gompertz","lognormal","uniform"),
       lty=1,col=c(1:4,7),xjust=1)
dev.off()
#=================================================================================
##lognormal
setEPS()
postscript(file="Fig18.eps")
theta_hat <- apply(matrix(unlist(sfLapply(1:tt, bupa3)),8,tt),1,mean)
par(mfrow=c(1,3))
#Density
plot(x=x,y=dlnorm(x,meanlog,sdlog),
     type="l",lty=1,col=1,ylab="Density",xlab="Time")
lines(x=x,y=dweibull(x,theta_hat[1],theta_hat[2]),type="l",lty=1,col=2)
lines(x=x,y=theta_hat[3]*exp(theta_hat[4]*x)*exp(-theta_hat[3]/theta_hat[4]*(exp(theta_hat[4]*x)-1)),type="l",lty=1,col=3)
lines(x=x,y=dlnorm(x,theta_hat[5],theta_hat[6]),type="l",lty=1,col=4)
lines(x=x,y=dunif(x,theta_hat[7],theta_hat[8]),type="l",lty=1,col=7)
legend("topright",legend=c("true","Weibull","Gompertz","lognormal","uniform"),
       lty=1,col=c(1:4,7),xjust=1)
#Survival
plot(x=x,y=plnorm(x,meanlog,sdlog,lower.tail = F),
     type="l",lty=1,col=1,ylab="Survival",xlab="Time")
lines(x=x,y=pweibull(x,theta_hat[1],theta_hat[2],lower.tail = F),type="l",lty=1,col=2)
lines(x=x,y=exp(-theta_hat[3]/theta_hat[4]*(exp(theta_hat[4]*x)-1)),type="l",lty=1,col=3)
lines(x=x,y=plnorm(x,theta_hat[5],theta_hat[6],lower.tail = F),type="l",lty=1,col=4)
lines(x=x,y=punif(x,theta_hat[7],theta_hat[8],lower.tail = F),type="l",lty=1,col=7)
legend("topright",legend=c("true","Weibull","Gompertz","lognormal","uniform"),
       lty=1,col=c(1:4,7),xjust=1)
#Hazard
plot(x=x,y=dlnorm(x,meanlog,sdlog)/plnorm(x,meanlog,sdlog,lower.tail = F),
     type="l",lty=1,col=1,ylab="Survival",xlab="Time")
lines(x=x,y=dweibull(x,theta_hat[1],theta_hat[2])/pweibull(x,theta_hat[1],theta_hat[2],lower.tail = F),type="l",lty=1,col=2)
lines(x=x,y=theta_hat[3]*exp(theta_hat[4]*x),type="l",lty=1,col=3)
lines(x=x,y=dlnorm(x,theta_hat[5],theta_hat[6])/plnorm(x,theta_hat[5],theta_hat[6],lower.tail = F),type="l",lty=1,col=4)
lines(x=x,y=dunif(x,theta_hat[7],theta_hat[8])/punif(x,theta_hat[7],theta_hat[8],lower.tail = F),type="l",lty=1,col=7)
legend("topright",legend=c("true","Weibull","Gompertz","lognormal","uniform"),
       lty=1,col=c(1:4,7),xjust=1)
dev.off()
#=================================================================================
##uniform
setEPS()
postscript(file="Fig19.eps")
theta_hat <- apply(matrix(unlist(sfLapply(1:tt, bupa4)),8,tt),1,mean)
sfStop()
par(mfrow=c(1,3))
#Density
plot(x=x,y=dunif(x,min_unif,max_unif),
     type="l",lty=1,col=1,ylab="Density",xlab="Time")
lines(x=x,y=dweibull(x,theta_hat[1],theta_hat[2]),type="l",lty=1,col=2)
lines(x=x,y=theta_hat[3]*exp(theta_hat[4]*x)*exp(-theta_hat[3]/theta_hat[4]*(exp(theta_hat[4]*x)-1)),type="l",lty=1,col=3)
lines(x=x,y=dlnorm(x,theta_hat[5],theta_hat[6]),type="l",lty=1,col=4)
lines(x=x,y=dunif(x,theta_hat[7],theta_hat[8]),type="l",lty=1,col=7)
legend("topright",legend=c("true","Weibull","Gompertz","lognormal","uniform"),
       lty=1,col=c(1:4,7),xjust=1)
#Survival
plot(x=x,y=punif(x,min_unif,max_unif,lower.tail = F),
     type="l",lty=1,col=1,ylab="Survival",xlab="Time")
lines(x=x,y=pweibull(x,theta_hat[1],theta_hat[2],lower.tail = F),type="l",lty=1,col=2)
lines(x=x,y=exp(-theta_hat[3]/theta_hat[4]*(exp(theta_hat[4]*x)-1)),type="l",lty=1,col=3)
lines(x=x,y=plnorm(x,theta_hat[5],theta_hat[6],lower.tail = F),type="l",lty=1,col=4)
lines(x=x,y=punif(x,theta_hat[7],theta_hat[8],lower.tail = F),type="l",lty=1,col=7)
legend("topright",legend=c("true","Weibull","Gompertz","lognormal","uniform"),
       lty=1,col=c(1:4,7),xjust=1)
#Hazard
plot(x=x,y=dunif(x,min_unif,max_unif)/punif(x,min_unif,max_unif,lower.tail = F),
     type="l",lty=1,col=1,ylab="Survival",xlab="Time")
lines(x=x,y=dweibull(x,theta_hat[1],theta_hat[2])/pweibull(x,theta_hat[1],theta_hat[2],lower.tail = F),type="l",lty=1,col=2)
lines(x=x,y=theta_hat[3]*exp(theta_hat[4]*x),type="l",lty=1,col=3)
lines(x=x,y=dlnorm(x,theta_hat[5],theta_hat[6])/plnorm(x,theta_hat[5],theta_hat[6],lower.tail = F),type="l",lty=1,col=4)
lines(x=x,y=dunif(x,theta_hat[7],theta_hat[8])/punif(x,theta_hat[7],theta_hat[8],lower.tail = F),type="l",lty=1,col=7)
legend("topright",legend=c("true","Weibull","Gompertz","lognormal","uniform"),
       lty=1,col=c(1:4,7),xjust=1)
dev.off()