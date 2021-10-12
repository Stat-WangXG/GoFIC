library(ReIns)
library(maxLik)
library(magrittr)
library(survival)
library(pracma)
library(doParallel)
library(snowfall)


source("f.R")
source("dlrb_n.R")
source("lamda_nbwB.R")
source("lamda_nbgB.R")
source("lamda_nblB.R")
source("Turnbull3.R")
source("gpz.R")
##=======================AIDS Clinical Trial Group study Data======================
B <- 5e2
aic <- bic <- hqc <- Lamda_n <- app_p <- rep(NA,3)
theta <- matrix(NA,3,2)
tsd <- matrix(NA,3,2)
#======================================================================
gof_ic <- function(data){
  n <- nrow(data)
  classb <- rep(1:B,each=n)
  data$censor <- data$L!=data$R
  alpha1 <- mean(data$censor==0)
  n1 <- round(n*alpha1)
  n2 <- round(n*(1-alpha1))
  m <- 1
  
  nu <- length(c(data$L,data$R[data$R<Inf]))
  U_data <- matrix(sort(c(data$L,data$R[data$R<Inf])),m,nu,byrow = T)
  non <- Turnbull(U_data,data$L,data$R,data$censor)
  tb2 <- non$surv
  phihat <- 1-tb2[nu]
  
  S1hat_time <- (phihat-1+tb2)/phihat
  dat2 <- data.frame(left=data$L,right=data$R,delta1=(data$L==0),delta2=(0<data$L&data$L<data$R&data$R<Inf),delta3=(data$R==Inf))
  
  ml_p1 <- (function(theta){f("weibull",dat2,theta,tb2)}) %>%
    maxNM(., start=c(.5,.5), constraints=list(ineqA= matrix(c(1,0,0,1),2,2), ineqB=rep(0,2)))
  aic[1] <- -2*ml_p1$maximum+4;  bic[1] <- -2*ml_p1$maximum+2*log(nrow(data));hqc[1] <- -2*ml_p1$maximum+2*log(log(nrow(data)))
  theta[1:2] <- theta_hat <- ml_p1$estimate
  S1_theta_time1 <- pweibull(q=U_data,shape=theta_hat[1],scale=theta_hat[2],lower.tail=FALSE)
  Lamda_n[1] <- sum((S1hat_time-S1_theta_time1)^2)/2
  Yb <- matrix(rbinom(n*B,1,phihat),n,B)
  Tb <- matrix(ifelse(Yb==1,rweibull(n*B,theta_hat[1],theta_hat[2]),Inf),n,B)
  sa_b <- sample(which(Yb==1),n1*B)
  Tb1 <- matrix(Tb[sa_b],n1,B)
  Tb2 <- matrix(Tb[-sa_b],n2,B)
  rm(Yb)
  Ub <- array(NA,c(m,n2,B))
  for(i in 1:m){
  Ub[i,,] <- sample(x=U_data[i,],size=n2*B,replace=TRUE)}
  datb_cb <- dlrbn(Tb1,Tb2,Ub,n,n1,m,B,classb)
  tb1_cb <- mapply(Turnbull3,datb_cb)
  lw <- lamda_nbwB(theta_hat,datb_cb,tb1_cb,B)
  tsd[1,] <- lw[[2]]
  app_p[1] <- sum(lw[[1]]>=Lamda_n[1])/B
  
  ml_p2 <- (function(theta){f("gompertz",dat2,theta,tb2)}) %>%
  maxNM(., start=c(0.1,0.1), constraints=list(ineqA= matrix(c(1,0,0,1),2,2), ineqB=rep(0,2)))
  aic[2] <- -2*ml_p2$maximum+4;bic[2] <- -2*ml_p2$maximum+2*log(nrow(data));hqc[2] <- -2*ml_p2$maximum+2*log(log(nrow(data)))
  theta[3:4] <- theta_hat <- ml_p2$estimate
  S1_theta_time2 <- pgpz(U_data,theta_hat[1],theta_hat[2],lower.tail = F)
  Lamda_n[2] <- sum((S1hat_time-S1_theta_time2)^2)/2
  Yb <- matrix(rbinom(n*B,1,phihat),n,B)
  Tb <- matrix(ifelse(Yb==1,rgpz(n*B,theta_hat[1],theta_hat[2]),Inf),n,B)
  sa_b <- sample(which(Yb==1),n1*B)
  Tb1 <- matrix(Tb[sa_b],n1,B)
  Tb2 <- matrix(Tb[-sa_b],n2,B)
  rm(Yb)
  Ub <- array(NA,c(m,n2,B))
  for(i in 1:m){
  Ub[i,,] <- sample(x=U_data[i,],size=n2*B,replace=TRUE)}
  datb_cb <- dlrbn(Tb1,Tb2,Ub,n,n1,m,B,classb)
  tb1_cb <- mapply(Turnbull3,datb_cb)
  lg <- lamda_nbgB(theta_hat,datb_cb,tb1_cb,B)
  tsd[2,] <- lg[[2]]
  app_p[2] <- sum(lg[[1]]>=Lamda_n[2])/B
  
  ml_p3 <- (function(theta){f("lognormal",dat2,theta,tb2)}) %>%
  maxNM(., start=c(1,2), constraints=list(ineqA= matrix(c(1,0,0,1),2,2), ineqB=rep(0,2)))
  aic[3] <- -2*ml_p3$maximum+4;bic[3] <- -2*ml_p3$maximum+2*log(nrow(data));hqc[3] <- -2*ml_p3$maximum+2*log(log(nrow(data)))
  theta[5:6] <- theta_hat <- ml_p3$estimate
  S1_theta_time3 <- plnorm(q=U_data,meanlog=theta_hat[1],sdlog=theta_hat[2],lower.tail=FALSE)
  Lamda_n[3] <- sum((S1hat_time-S1_theta_time3)^2)/2
  Yb <- matrix(rbinom(n*B,1,phihat),n,B)
  Tb <- matrix(ifelse(Yb==1,rlnorm(n*B,theta_hat[1],theta_hat[2]),Inf),n,B)
  sa_b <- sample(which(Yb==1),n1*B)
  Tb1 <- matrix(Tb[sa_b],n1,B)
  Tb2 <- matrix(Tb[-sa_b],n2,B)
  rm(Yb)
  Ub <- array(NA,c(m,n2,B))
  for(i in 1:m){
    Ub[i,,] <- sample(x=U_data[i,],size=n2*B,replace=TRUE)}
  datb_cb <- dlrbn(Tb1,Tb2,Ub,n,n1,m,B,classb)
  tb1_cb <- mapply(Turnbull3,datb_cb)
  ll <- lamda_nblB(theta_hat,datb_cb,tb1_cb,B)
  tsd[3,] <- ll[[2]]
  app_p[3] <- sum(ll[[1]]>=Lamda_n[3])/B
  
return(list(p_value=app_p,theta=theta,s.d._of_theta=tsd,AIC=aic))}
