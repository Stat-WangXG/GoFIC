myfun3 <- function(aaa){
  tic("total")
  ##=============================caculation========================
  tic("preparation")
  T3 <- T1_func("weibull")
  T4 <- c(T3, rep(Inf, round(n4))) #length:n
  
  sa <- sample(n3,n1)
  T1 <- T3[sa]
  T2 <- T4[-sa]
  U <- matrix(NA,m,n2)
  U[1,] <- rweibull(n2,1,0.5)
  for(i in 2:m){
    U[i,] <- max(U[i-1,]) + rweibull(n2,1,0.5)}
  
  rm(T3)
  dat <- dlr(T1,T2,U)[[1]]
  x <- sort(c(dat$left[dat$left>0],dat$right[dat$right<Inf]))
  tb1 <- Turnbull(x,dat$left,dat$right,rowSums(dat[,3:5]))$surv
  #========================================lognormal====================================
  #============Parametric logLikelihood(get theta hat)==============
  ml_p <- (function(theta){f("lognormal",dat,theta,tb1)}) %>%
    maxLik(.,start=c(0.2,1),method = "NM",constraints = list(ineqA=matrix(c(1,0,0,1),2,2), ineqB=c(-1e-6,-1e-6)))
  
  theta_hat <- ml_p$estimate
  toc()
  #====================Leveraged Bootstrap=======================
  tic("le")
  rej <- leverage_lo(x,tb1,theta_hat)
  toc()
  tic("le1")
  rejx <- leverage_lo1(x,tb1,theta_hat)
  toc()
  #==============================================================
  #========================return============================
  return(list(rej,rejx))
  toc()}