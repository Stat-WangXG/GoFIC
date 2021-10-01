we_150 <- read.csv("C:/Users/Administrator/Documents/Statistics/待议论文/cure/GOF/GOF for Turnbull/tttttry/150/alpha0.1phi0.5(2021)/we.csv",header = F)
go_150 <- read.csv("C:/Users/Administrator/Documents/Statistics/待议论文/cure/GOF/GOF for Turnbull/tttttry/150/alpha0.1phi0.5(2021)/go.csv",header = F)
lo_150 <- read.csv("C:/Users/Administrator/Documents/Statistics/待议论文/cure/GOF/GOF for Turnbull/tttttry/150/alpha0.1phi0.5(2021)/lo.csv",header = F)
un_150 <- read.csv("C:/Users/Administrator/Documents/Statistics/待议论文/cure/GOF/GOF for Turnbull/tttttry/150/alpha0.1phi0.5(2021)/un.csv",header = F)

we_500 <- read.csv("C:/Users/Administrator/Documents/Statistics/待议论文/cure/GOF/GOF for Turnbull/tttttry/500/alpha0.1phi0.5(2021)/we.csv",header = F)
go_500 <- read.csv("C:/Users/Administrator/Documents/Statistics/待议论文/cure/GOF/GOF for Turnbull/tttttry/500/alpha0.1phi0.5(2021)/go.csv",header = F)
lo_500 <- read.csv("C:/Users/Administrator/Documents/Statistics/待议论文/cure/GOF/GOF for Turnbull/tttttry/500/alpha0.1phi0.5(2021)/lo.csv",header = F)
un_500 <- read.csv("C:/Users/Administrator/Documents/Statistics/待议论文/cure/GOF/GOF for Turnbull/tttttry/500/alpha0.1phi0.5(2021)/un.csv",header = F)

write.table(round(matrix(c(apply(we_150[,16:24],2,mean),apply(we_150[,16:24],2,sd),
                           apply(we_500[,16:24],2,mean),apply(we_500[,16:24],2,sd),
                           apply(go_150[,16:24],2,mean),apply(go_150[,16:24],2,sd),
                           apply(go_500[,16:24],2,mean),apply(go_500[,16:24],2,sd),
                           apply(lo_150[,16:24],2,mean),apply(lo_150[,16:24],2,sd),
                           apply(lo_500[,16:24],2,mean),apply(lo_500[,16:24],2,sd),
                           apply(un_150[,16:24],2,mean),apply(un_150[,16:24],2,sd),
                           apply(un_500[,16:24],2,mean),apply(un_500[,16:24],2,sd)
),9,16),3),file="附录Table5.csv",
row.names = F,col.names = c("we_150m","we_150s","we_500m","we_500s","go_150m","go_150s","go_500m","go_500s","lo_150m","lo_150s","lo_500m","lo_500s","un_150m","un_150s","un_500m","un_500s"),sep = ",")
