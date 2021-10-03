myfunc4 <- function(aaa){
  
  ##=============================caculation========================
  tic("total")
  T3 <- T1_func("we+20lo")
  T4 <- c(T3, rep(Inf, round(n4)))
  sa <- sample(n3,n1)
  T1 <- T3[sa]
  T2 <- T4[-sa]
  
  U <- matrix(NA,m,n2)
  U[1,] <- rweibull(n2,1,0.5)
  for(i in 2:m){
    U[i,] <- max(U[i-1,]) + rweibull(n2,1,0.5)}
  
  rm(T3,T4)
  tic("NPMLE")
  dlr_ <- dlr(T1,T2,U)
  dat <- dlr_[[1]]
  ic <- dlr_[[2]]
  rc <- dlr_[[3]]
  
  tb1 <- Turnbull3(dat)
  #tb2 <- Turnbull(c(dat$left[dat$delta3,],dat$right[dat$delta1,]),dat$left,dat$right,censored = rep(1,sum(dat$delta2==1))$surv
  phihat <- 1-min(tb1)
  
  S1hat_time1 <- (phihat-1+tb1)/phihat
  #Turnbull估计得到的基于所观察到数据的易感子群体
  S1_time1 <- pweibull(q=c(dat$left[dat$left>0],dat$right[dat$right<Inf]),shape=shape_weib,scale=scale_weib,lower.tail=FALSE) 
  #真正的易感子群体
  MSE_np <- sum((S1_time1-S1hat_time1)^2)/2#+sum((S1_time2-S1hat_time2)^2)/2
  toc()
  #====================================Weibull====================================
  
  tic("PMLE for Weibull")
  ml_p <- (function(theta){f("weibull",dat,theta,tb1)}) %>%
    maxLik(.,start=c(1,1),method = "NM",constraints = list(ineqA=matrix(c(1,0,0,1),2,2), ineqB=c(-1e-6,-1e-6)))
  aic <- -2*ml_p$maximum+4
  theta_hat <- ml_p$estimate
  #=================Test statitic(get Lambda_n)====================
  S1_theta_time1 <- pweibull(q=c(dat$left[dat$left>0],dat$right[dat$right<Inf]),shape=theta_hat[1],scale=theta_hat[2],lower.tail=FALSE)
  #S1_theta_time2 <- pweibull(q=c(dat$left[dat$delta3],dat$right[dat$delta1]),shape=theta_hat[1],scale=theta_hat[2],lower.tail=FALSE)
  MSE_p <- sum((S1_time1-S1_theta_time1)^2)/2#+sum((S1_time2-S1_theta_time2)^2)/2
  Lamda_n <- sum((S1hat_time1-S1_theta_time1)^2)/2#+sum((S1hat_time2-S1_theta_time2)^2)/2
  toc()
  #=================Bootstrap(get Lambda_n^b)====================
  tic("bootstrap data generation(W)")
  Yb <- matrix(rbinom(n*B,1,phihat),n,B)
  Tb <- matrix(ifelse(Yb==1,rweibull(n*B,theta_hat[1],theta_hat[2]),Inf),n,B)
  
  #==========================sa,U,data,tb1========================
  sa_b <- sample(which(Yb==1),n1*B)
  Tb1 <- matrix(Tb[sa_b],n1,B)
  Tb2 <- matrix(Tb[-sa_b],n2,B)
  rm(Yb)
  
  Ub <- array(NA,c(m,n2,B))
  for(i in 1:m){
    Ub[i,,] <- sample(x=U[i,],size=n2*B,replace=TRUE)
  }
  #for(i in 1:(m-1)){
  #Ub[i,,] <- sample(x=U[i,],size=n2*B,replace=TRUE)
  #}
  #for(b in 1:B){
  #Ub[m,,b] <- max(c(Tb1[,b],Tb2[,b][Tb2[,b]<Inf],U[m-1,])) + 1
  #}
  toc()
  tic("dlrb(W)")
  datb_cb <- dlrb(Tb1,Tb2,Ub)
  toc()
  tic("Turnbull3(W)")
  tb1_cb <- mapply(Turnbull3,datb_cb)
  toc()
  #===========================approximate p value =================================
  tic("nbw")
  app_p <- sum(lamda_nbw(theta_hat,datb_cb,tb1_cb)[[1]]>=Lamda_n)/B
  toc()
  return(list(app_p,MSE_np,MSE_p,ic,rc,aic))
  toc()}