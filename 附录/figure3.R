we <- read.csv("C:/Users/Administrator/Documents/Statistics/待议论文/conditional cure延伸/cure/GOF/GOF for Turnbull/tttttry/500co/alpha0.1phi0.8/we.csv",header = F)
go <- read.csv("C:/Users/Administrator/Documents/Statistics/待议论文/conditional cure延伸/cure/GOF/GOF for Turnbull/tttttry/500co/alpha0.1phi0.8/go.csv",header = F)
lo <- read.csv("C:/Users/Administrator/Documents/Statistics/待议论文/conditional cure延伸/cure/GOF/GOF for Turnbull/tttttry/500co/alpha0.1phi0.8/lo.csv",header = F)
w5g <- read.csv("C:/Users/Administrator/Documents/Statistics/待议论文/conditional cure延伸/cure/GOF/GOF for Turnbull/tttttry/500co/alpha0.1phi0.8/we+5go.csv",header = F)
w10g <- read.csv("C:/Users/Administrator/Documents/Statistics/待议论文/conditional cure延伸/cure/GOF/GOF for Turnbull/tttttry/500co/alpha0.1phi0.8/we+10go.csv",header = F)
w15g <- read.csv("C:/Users/Administrator/Documents/Statistics/待议论文/conditional cure延伸/cure/GOF/GOF for Turnbull/tttttry/500co/alpha0.1phi0.8/we+15go.csv",header = F)
w20g <- read.csv("C:/Users/Administrator/Documents/Statistics/待议论文/conditional cure延伸/cure/GOF/GOF for Turnbull/tttttry/500co/alpha0.1phi0.8/we+20go.csv",header = F)
w5l <- read.csv("C:/Users/Administrator/Documents/Statistics/待议论文/conditional cure延伸/cure/GOF/GOF for Turnbull/tttttry/500co/alpha0.1phi0.8/we+5lo.csv",header = F)
w10l <- read.csv("C:/Users/Administrator/Documents/Statistics/待议论文/conditional cure延伸/cure/GOF/GOF for Turnbull/tttttry/500co/alpha0.1phi0.8/we+10lo.csv",header = F)
w15l <- read.csv("C:/Users/Administrator/Documents/Statistics/待议论文/conditional cure延伸/cure/GOF/GOF for Turnbull/tttttry/500co/alpha0.1phi0.8/we+15lo.csv",header = F)
w20l <- read.csv("C:/Users/Administrator/Documents/Statistics/待议论文/conditional cure延伸/cure/GOF/GOF for Turnbull/tttttry/500co/alpha0.1phi0.8/we+20lo.csv",header = F)

#jpeg(file = "fig1.jpg",width=1200,height=1200)
setEPS()
postscript(file="Fig3.eps")
par(mfrow=c(2,2))
plot(c(seq(0,20,5),sqrt(500)),c(mean(we[1,]<.05),mean(w5g[1,]<.05),mean(w10g[1,]<.05),mean(w15g[1,]<.05),mean(w20g[1,]<.05),mean(go[1,]<.05)),
     type = "o",pch=1,xlim=c(0,25),ylim = c(0,1),xlab = "c",ylab = "")
lines(c(seq(0,20,5),sqrt(500)),c(mean(we[1,]<0.1),mean(w5g[1,]<0.1),mean(w10g[1,]<0.1),mean(w15g[1,]<0.1),mean(w20g[1,]<0.1),mean(go[1,]<0.1)),type = "o",pch=2)
lines(c(seq(0,20,5),sqrt(500)),c(1,apply(w5g[6,],1,mean)/apply(we[6,],1,mean),apply(w10g[6,],1,mean)/apply(we[6,],1,mean),apply(w15g[6,],1,mean)/apply(we[6,],1,mean),
                                 apply(w20g[6,],1,mean)/apply(we[6,],1,mean),apply(go[6,],1,mean)/apply(we[6,],1,mean)),type = "o",pch=4)

plot(c(seq(0,20,5),sqrt(500)),c(mean(we[1,]<.05),mean(w5l[1,]<.05),mean(w10l[1,]<.05),mean(w15l[1,]<.05),mean(w20l[1,]<.05),mean(lo[1,]<.05)),
     type = "o",pch=1,xlim=c(0,25),ylim = c(0,1),xlab = "c",ylab = "")
lines(c(seq(0,20,5),sqrt(500)),c(mean(we[1,]<0.1),mean(w5l[1,]<0.1),mean(w10l[1,]<0.1),mean(w15l[1,]<0.1),mean(w20l[1,]<0.1),mean(lo[1,]<0.1)),type = "o",pch=2)
lines(c(seq(0,20,5),sqrt(500)),c(1,apply(w5l[6,],1,mean)/apply(we[6,],1,mean),apply(w10l[6,],1,mean)/apply(we[6,],1,mean),apply(w15l[6,],1,mean)/apply(we[6,],1,mean),
                                 apply(w20l[6,],1,mean)/apply(we[6,],1,mean),apply(lo[6,],1,mean)/apply(we[6,],1,mean)),type = "o",pch=4)

we_ <- read.csv("C:/Users/Administrator/Documents/Statistics/待议论文/conditional cure延伸/cure/GOF/GOF for Turnbull/tttttry/500co/alpha0.3phi0.5/we.csv",header = F)
go_ <- read.csv("C:/Users/Administrator/Documents/Statistics/待议论文/conditional cure延伸/cure/GOF/GOF for Turnbull/tttttry/500co/alpha0.3phi0.5/go.csv",header = F)
lo_ <- read.csv("C:/Users/Administrator/Documents/Statistics/待议论文/conditional cure延伸/cure/GOF/GOF for Turnbull/tttttry/500co/alpha0.3phi0.5/lo.csv",header = F)
w5g_ <- read.csv("C:/Users/Administrator/Documents/Statistics/待议论文/conditional cure延伸/cure/GOF/GOF for Turnbull/tttttry/500co/alpha0.3phi0.5/we+5go.csv",header = F)
w10g_ <- read.csv("C:/Users/Administrator/Documents/Statistics/待议论文/conditional cure延伸/cure/GOF/GOF for Turnbull/tttttry/500co/alpha0.3phi0.5/we+10go.csv",header = F)
w15g_ <- read.csv("C:/Users/Administrator/Documents/Statistics/待议论文/conditional cure延伸/cure/GOF/GOF for Turnbull/tttttry/500co/alpha0.3phi0.5/we+15go.csv",header = F)
w20g_ <- read.csv("C:/Users/Administrator/Documents/Statistics/待议论文/conditional cure延伸/cure/GOF/GOF for Turnbull/tttttry/500co/alpha0.3phi0.5/we+20go.csv",header = F)
w5l_ <- read.csv("C:/Users/Administrator/Documents/Statistics/待议论文/conditional cure延伸/cure/GOF/GOF for Turnbull/tttttry/500co/alpha0.3phi0.5/we+5lo.csv",header = F)
w10l_ <- read.csv("C:/Users/Administrator/Documents/Statistics/待议论文/conditional cure延伸/cure/GOF/GOF for Turnbull/tttttry/500co/alpha0.3phi0.5/we+10lo.csv",header = F)
w15l_ <- read.csv("C:/Users/Administrator/Documents/Statistics/待议论文/conditional cure延伸/cure/GOF/GOF for Turnbull/tttttry/500co/alpha0.3phi0.5/we+15lo.csv",header = F)
w20l_ <- read.csv("C:/Users/Administrator/Documents/Statistics/待议论文/conditional cure延伸/cure/GOF/GOF for Turnbull/tttttry/500co/alpha0.3phi0.5/we+20lo.csv",header = F)

plot(c(seq(0,20,5),sqrt(500)),c(mean(we_[1,]<.05),mean(w5g_[1,]<.05),mean(w10g_[1,]<.05),mean(w15g_[1,]<.05),mean(w20g_[1,]<.05),mean(go_[1,]<.05)),
     type = "o",pch=1,xlim=c(0,25),ylim = c(0,1),xlab = "c",ylab = "")
lines(c(seq(0,20,5),sqrt(500)),c(mean(we_[1,]<0.1),mean(w5g_[1,]<0.1),mean(w10g_[1,]<0.1),mean(w15g_[1,]<0.1),mean(w20g_[1,]<0.1),mean(go_[1,]<0.1)),type = "o",pch=2)
lines(c(seq(0,20,5),sqrt(500)),c(1,apply(w5g_[12,],1,mean)/apply(we_[6,],1,mean),apply(w10g_[12,],1,mean)/apply(we_[6,],1,mean),apply(w15g_[12,],1,mean)/apply(we_[6,],1,mean),
                                 apply(w20g_[12,],1,mean)/apply(we_[6,],1,mean),apply(go_[6,],1,mean)/apply(we_[6,],1,mean)),type = "o",pch=4)

plot(c(seq(0,20,5),sqrt(500)),c(mean(we_[1,]<.05),mean(w5l_[1,]<.05),mean(w10l_[1,]<.05),mean(w15l_[1,]<.05),mean(w20l_[1,]<.05),mean(lo_[1,]<.05)),
     type = "o",pch=1,xlim=c(0,25),ylim = c(0,1),xlab = "c",ylab = "")
lines(c(seq(0,20,5),sqrt(500)),c(mean(we_[1,]<0.1),mean(w5l_[1,]<0.1),mean(w10l_[1,]<0.1),mean(w15l_[1,]<0.1),mean(w20l_[1,]<0.1),mean(lo_[1,]<0.1)),type = "o",pch=2)
lines(c(seq(0,20,5),sqrt(500)),c(1,apply(w5l_[12,],1,mean)/apply(we_[6,],1,mean),apply(w10l_[12,],1,mean)/apply(we_[6,],1,mean),apply(w15l_[12,],1,mean)/apply(we_[6,],1,mean),
                                 apply(w20l_[12,],1,mean)/apply(we_[6,],1,mean),apply(lo_[6,],1,mean)/apply(we_[6,],1,mean)),type = "o",pch=4)
dev.off()
