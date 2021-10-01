we <- read.csv("C:/Users/Administrator/Documents/Statistics/待议论文/cure/GOF/GOF for Turnbull/tttttry/500/alpha0.3phi0.8(2021)/we.csv",header = F)
go <- read.csv("C:/Users/Administrator/Documents/Statistics/待议论文/cure/GOF/GOF for Turnbull/tttttry/500/alpha0.3phi0.8(2021)/go.csv",header = F)
lo <- read.csv("C:/Users/Administrator/Documents/Statistics/待议论文/cure/GOF/GOF for Turnbull/tttttry/500/alpha0.3phi0.8(2021)/lo.csv",header = F)
w5g <- read.csv("C:/Users/Administrator/Documents/Statistics/待议论文/cure/GOF/GOF for Turnbull/tttttry/500co/alpha0.3phi0.8(2021)/we+5go.csv",header = F)
w10g <- read.csv("C:/Users/Administrator/Documents/Statistics/待议论文/cure/GOF/GOF for Turnbull/tttttry/500co/alpha0.3phi0.8(2021)/we+10go.csv",header = F)
w15g <- read.csv("C:/Users/Administrator/Documents/Statistics/待议论文/cure/GOF/GOF for Turnbull/tttttry/500co/alpha0.3phi0.8(2021)/we+15go.csv",header = F)
w20g <- read.csv("C:/Users/Administrator/Documents/Statistics/待议论文/cure/GOF/GOF for Turnbull/tttttry/500co/alpha0.3phi0.8(2021)/we+20go.csv",header = F)
w5l <- read.csv("C:/Users/Administrator/Documents/Statistics/待议论文/cure/GOF/GOF for Turnbull/tttttry/500co/alpha0.3phi0.8(2021)/we+5lo.csv",header = F)
w10l <- read.csv("C:/Users/Administrator/Documents/Statistics/待议论文/cure/GOF/GOF for Turnbull/tttttry/500co/alpha0.3phi0.8(2021)/we+10lo.csv",header = F)
w15l <- read.csv("C:/Users/Administrator/Documents/Statistics/待议论文/cure/GOF/GOF for Turnbull/tttttry/500co/alpha0.3phi0.8(2021)/we+15lo.csv",header = F)
w20l <- read.csv("C:/Users/Administrator/Documents/Statistics/待议论文/cure/GOF/GOF for Turnbull/tttttry/500co/alpha0.3phi0.8(2021)/we+20lo.csv",header = F)


write.table(matrix(c(apply(we[,10:11],2,mean),apply(we[,5:6],2,mean),mean(we[,1]<.05),mean(we[,1]<.1),mean(we[,12]),
                     apply(w5g[,4:5],2,mean),apply(w5g[,2:3],2,mean),mean(w5g[,1]<.05),mean(w5g[,1]<.1),mean(w5g[,6]),
                     apply(w10g[,4:5],2,mean),apply(w10g[,2:3],2,mean),mean(w10g[,1]<.05),mean(w10g[,1]<.1),mean(w10g[,6]),
                     apply(w15g[,4:5],2,mean),apply(w15g[,2:3],2,mean),mean(w15g[,1]<.05),mean(w15g[,1]<.1),mean(w15g[,6]),
                     apply(w20g[,4:5],2,mean),apply(w20g[,2:3],2,mean),mean(w20g[,1]<.05),mean(w20g[,1]<.1),mean(w20g[,6]),
                     apply(go[,10:11],2,mean),apply(go[,5:6],2,mean),mean(go[,1]<.05),mean(go[,1]<.1),mean(go[,12]),
                     apply(we[,10:11],2,mean),apply(we[,5:6],2,mean),mean(we[,1]<.05),mean(we[,1]<.1),mean(we[,12]),
                     apply(w5l[,4:5],2,mean),apply(w5l[,2:3],2,mean),mean(w5l[,1]<.05),mean(w5l[,1]<.1),mean(w5l[,6]),
                     apply(w10l[,4:5],2,mean),apply(w10l[,2:3],2,mean),mean(w10l[,1]<.05),mean(w10l[,1]<.1),mean(w10l[,6]),
                     apply(w15l[,4:5],2,mean),apply(w15l[,2:3],2,mean),mean(w15l[,1]<.05),mean(w15l[,1]<.1),mean(w15l[,6]),
                     apply(w20l[,4:5],2,mean),apply(w20l[,2:3],2,mean),mean(w20l[,1]<.05),mean(w20l[,1]<.1),mean(w20l[,6]),
                     apply(lo[,10:11],2,mean),apply(lo[,5:6],2,mean),mean(lo[,1]<.05),mean(lo[,1]<.1),mean(lo[,12])
                     ),12,7,byrow = T),file="附录Table2.csv",row.names = F,col.names = c("ic","rc","RMSEnp","RMSEp","r5","r10","AIC"),sep = ",")

list(apply(lo[,4:5],2,mean),apply(lo[,2:3],2,mean),mean(lo[,1]<.05),mean(lo[,1]<.1),apply(lo[,6]))
list(apply(go[,4:5],2,mean),apply(go[,2:3],2,mean),mean(go[,1]<.05),mean(go[,1]<.1),apply(go[,6]))
list(apply(we[,4:5],2,mean),apply(we[,2:3],2,mean),mean(we[,1]<.05),mean(we[,1]<.1),apply(we[,6]))