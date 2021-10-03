myfun <- function(aaa){
  
  ##=============================caculation========================
  tic("total")
  aic <- rep(NA,4)
  app_p <- rep(NA,4)
  MSE_p <- rep(NA,4)
  theta <- rep(NA,8)
  
  T3 <- T1_func("weibull")
  T4 <- c(T3, rep(Inf, round(n4))) #length:n
  
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
  
  #############             #############         #############         #############
  #############             #############         #############         #############
  #############             #############         #############         #############
  #############             #############         #############         #############
  S1hat_time1 <- (phihat-1+tb1)/phihat
  #S1hat_time2 <- (phihat-1+tb2)/phihat
  #Turnbull估计得到的基于所观察到数据的易感子群体
  S1_time1 <- pweibull(q=c(dat$left[dat$left>0],dat$right[dat$right<Inf]),shape=shape_weib,scale=scale_weib,lower.tail=FALSE) 
  #S1_time2 <- pweibull(q=c(dat$left[dat$delta3],dat$right[dat$delta1]),shape=shape_weib,scale=scale_weib,lower.tail=FALSE) 
  
  #真正的易感子群体
  MSE_np <- sum((S1_time1-S1hat_time1)^2)/2#+sum((S1_time2-S1hat_time2)^2)/2
  toc()
  #====================================Weibull====================================
  tic("PMLE for Weibull")
  ml_p <- (function(theta){f("weibull",dat,theta,tb1)}) %>%
    maxLik(.,start=c(1,1),method = "NM",constraints = list(ineqA=matrix(c(1,0,0,1),2,2), ineqB=c(-1e-6,-1e-6)))
  aic[1] <- -2*ml_p$maximum+4
  theta_hat <- ml_p$estimate
  theta[1:2] <- theta_hat
  #=================Test statitic(get Lambda_n)====================
  S1_theta_time1 <- pweibull(q=c(dat$left[dat$left>0],dat$right[dat$right<Inf]),shape=theta_hat[1],scale=theta_hat[2],lower.tail=FALSE)
  #S1_theta_time2 <- pweibull(q=c(dat$left[dat$delta3],dat$right[dat$delta1]),shape=theta_hat[1],scale=theta_hat[2],lower.tail=FALSE)
  MSE_p[1] <- sum((S1_time1-S1_theta_time1)^2)/2#+sum((S1_time2-S1_theta_time2)^2)/2
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
  lbw <- lamda_nbw(theta_hat,datb_cb,tb1_cb)
  toc()
  app_p[1] <- lbw[[1]] %>% (function(x){sum(x>=Lamda_n)/B})
  
  ##########===============Uniform===================================
  #============Parametric logLikelihood(get theta hat)==============
  tic("PMLE for uniform")
  ml_p <- (function(theta){f("uniform",dat,theta,tb1)}) %>%
    maxLik(.,start=c(0,max(dat$right[dat$right<Inf])),method = "NM",constraints = list(ineqA=matrix(c(1,-1,0,1),2,2), ineqB=c(0.3e-323,-1e-6)))
  
  aic[4] <- -2*ml_p$maximum+4
  theta_hat <- ml_p$estimate
  theta[7:8] <- theta_hat
  
  #=================Test statitic(get Lambda_n)====================
  
  S1_theta_time4 <- punif(q=c(dat$left[dat$left>0],dat$right[dat$right<Inf]),min=theta_hat[1],max=theta_hat[2],lower.tail=FALSE)
  MSE_p[4] <- sum((S1_time1-S1_theta_time4)^2)/2
  Lamda_n <- sum((S1hat_time1-S1_theta_time4)^2)/2
  toc()
  #=============================Bootstrap(get Lambda_n^b)================================
  tic("bootstrap data generation(u)")
  Yb <- matrix(rbinom(n*B,1,phihat),n,B)
  Tb <- matrix(ifelse(Yb==1,runif(n*B,max(c(theta_hat[1],0)),theta_hat[2]),Inf),n,B)
  
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
  tic("dlrb(u)")
  datb_cb <- dlrb(Tb1,Tb2,Ub)
  toc()
  tic("Turnbull3(u)")
  tb1_cb <- mapply(Turnbull3,datb_cb)
  toc()
  #===========================approximate p value =================================
  tic("nbu")
  lbu <- lamda_nbu(theta_hat,datb_cb,tb1_cb)
  toc()
  app_p[4] <- lbu[[1]] %>% (function(x){sum(x>=Lamda_n)/B})
  
  
  
  #========================================lognormal====================================
  #============Parametric logLikelihood(get theta hat)==============
  tic("PMLE for lognormal")
  ml_p <- (function(theta){f("lognormal",dat,theta,tb1)}) %>%
    maxLik(.,start=c(1e-5,1),method = "NM",constraints = list(ineqA=matrix(c(1,0,0,1),2,2), ineqB=c(-1e-6,-1e-6)))
  
  aic[3] <- -2*ml_p$maximum+4
  theta_hat <- ml_p$estimate
  theta[5:6] <- theta_hat
  
  #=================Test statitic(get Lambda_n)====================
  S1_theta_time3 <- plnorm(q=c(dat$left[dat$left>0],dat$right[dat$right<Inf]),meanlog=theta_hat[1],sdlog=theta_hat[2],lower.tail=FALSE)
  MSE_p[3] <- sum((S1_time1-S1_theta_time3)^2)/2
  Lamda_n <- sum((S1hat_time1-S1_theta_time3)^2)/2
  toc()
  #=================Bootstrap(get Lambda_n^b)====================
  tic("bootstrap data generation(l)")
  Yb <- matrix(rbinom(n*B,1,phihat),n,B)
  Tb <- matrix(ifelse(Yb==1,rlnorm(n*B,theta_hat[1],theta_hat[2]),Inf),n,B)
  
  #==========================sa,U,data,tb1========================
  sa_b <- sample(which(Yb==1),n1*B)
  Tb1 <- matrix(Tb[sa_b],n1,B)
  Tb2 <- matrix(Tb[-sa_b],n2,B)
  rm(Yb)
  
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
  tic("dlrb(l)")
  datb_cb <- dlrb(Tb1,Tb2,Ub)
  toc()
  tic("Turnbull3(l)")
  tb1_cb <- mapply(Turnbull3,datb_cb)
  toc()
  #===========================approximate p value =================================
  tic("nbl")
  lbl <- lamda_nbl(theta_hat,datb_cb,tb1_cb)
  toc()
  app_p[3] <- lbl[[1]] %>% (function(x){sum(x>=Lamda_n)/B})
  
  
  
  #====================================Gompertz====================================
  tic("PMLE for Gompertz")
  ml_p <- (function(theta){f("gompertz",dat,theta,tb1)}) %>%
    maxLik(.,start=c(0.1,0.5),method = "NM",constraints = list(ineqA=matrix(c(1,0,0,1),2,2), ineqB=c(-1e-6,-1e-6)))
  
  aic[2] <- -2*ml_p$maximum+4
  theta_hat <- ml_p$estimate
  theta[3:4] <- theta_hat
  
  #=================Test statitic(get Lambda_n)====================
  S1_theta_time2 <- exp(-theta_hat[1]/theta_hat[2]*(exp(theta_hat[2]*c(dat$left[dat$left>0],dat$right[dat$right<Inf]))-1))
  MSE_p[2] <- sum((S1_time1-S1_theta_time2)^2)/2
  Lamda_n <- sum((S1hat_time1-S1_theta_time2)^2)/2
  toc()
  #=================Bootstrap(get Lambda_n^b)====================
  tic("bootstrap data generation(G)")
  Yb <- matrix(rbinom(n*B,1,phihat),n,B)
  Tb <- matrix(ifelse(Yb==1,rgpz(n*B,theta_hat[1],theta_hat[2]),Inf),n,B)
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
  tic("dlrb(G)")
  datb_cb <- dlrb(Tb1,Tb2,Ub)
  toc()
  tic("Turnbull3(G)")
  tb1_cb <- mapply(Turnbull3,datb_cb)
  toc()
  #===========================approximate p value =================================
  tic("nbg")
  lbg <- lamda_nbg(theta_hat,datb_cb,tb1_cb)
  app_p[2] <- lbg[[1]] %>% (function(x){sum(x>=Lamda_n)/B})
  toc()
  #========================return============================
  toc()
  return(list(app_p,MSE_np,MSE_p,ic,rc,aic,theta,phihat))
  
}