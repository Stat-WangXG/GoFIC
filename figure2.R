source("ini_value.R")
setEPS()
postscript(file="C:/Users/Administrator/Documents/Statistics/´ıÒéÂÛÎÄ/cure/GOF/GOF for Turnbull/tex/Fig2(1).eps")
par(mfrow=c(1,2))
x <- seq(0,8,1e-3)
plot(x=x,y=pweibull(x,shape_weib,scale_weib,lower.tail = F),
     type="l",lty=1,lwd=1,ylab="Survival",xlab="Time")
lines(x=x,y=pweibull(x,shape_weib,scale_weib,lower.tail=F) - sqrt(5)/10*(pweibull(x,shape_weib,scale_weib,lower.tail=F) - exp(-(labda*(exp(gama*x) - 1))/gama)),
      type="l",lty=2,lwd=1)
lines(x=x,y=pweibull(x,shape_weib,scale_weib,lower.tail=F) - 1/sqrt(5)*(pweibull(x,shape_weib,scale_weib,lower.tail=F) - exp(-(labda*(exp(gama*x) - 1))/gama)),
      type="l",lty=3,lwd=1)
lines(x=x,y=pweibull(x,shape_weib,scale_weib,lower.tail=F) - 3*5^(1/2)/10*(pweibull(x,shape_weib,scale_weib,lower.tail=F) - exp(-(labda*(exp(gama*x) - 1))/gama)),
      type="l",lty=4,lwd=1)
lines(x=x,y=pweibull(x,shape_weib,scale_weib,lower.tail=F) - 4*5^(1/2)/10*(pweibull(x,shape_weib,scale_weib,lower.tail=F) - exp(-(labda*(exp(gama*x) - 1))/gama)),
      type="l",lty=5,lwd=1)
lines(x=x,y=exp(-labda*(exp(gama*x) - 1)/gama),
      type="l",lty=6,lwd=1)
legend("topright",cex=0.8,legend=c("c=0","c=5","c=10","c=15","c=20","c=sqrt(500)"),
       lty=1:5)
plot(x=x,y=pweibull(x,shape_weib,scale_weib,lower.tail = F),
     type="l",lty=1,lwd=1,ylab="Survival",xlab="Time")
lines(x=x,y=pweibull(x,shape_weib,scale_weib,lower.tail=F) - sqrt(5)/10*(pweibull(x,shape_weib,scale_weib,lower.tail=F) - plnorm(x,meanlog,sdlog,lower.tail=FALSE)),
      type="l",lty=2,lwd=1)
lines(x=x,y=pweibull(x,shape_weib,scale_weib,lower.tail=F) - 1/sqrt(5)*(pweibull(x,shape_weib,scale_weib,lower.tail=F) - plnorm(x,meanlog,sdlog,lower.tail=FALSE)),
      type="l",lty=3,lwd=1)
lines(x=x,y=pweibull(x,shape_weib,scale_weib,lower.tail=F) - 3*5^(1/2)/10*(pweibull(x,shape_weib,scale_weib,lower.tail=F) - plnorm(x,meanlog,sdlog,lower.tail=FALSE)),
      type="l",lty=4,lwd=1)
lines(x=x,y=pweibull(x,shape_weib,scale_weib,lower.tail=F) - 4*5^(1/2)/10*(pweibull(x,shape_weib,scale_weib,lower.tail=F) - plnorm(x,meanlog,sdlog,lower.tail=FALSE)),
      type="l",lty=5,lwd=1)
lines(x=x,y=exp(-labda*(exp(gama*x) - 1)/gama),
      type="l",lty=6,lwd=1)
legend("topright",cex=0.8,legend=c("c=0","c=5","c=10","c=15","c=20","c=sqrt(500)"),
       lty=1:5,xjust=1)
dev.off()