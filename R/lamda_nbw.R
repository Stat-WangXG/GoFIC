lamda_nbw <- function(theta_hat,datb_cb,tb1_cb){
  theta_hatb <- matrix(NA,2,B)
  Lamda_nb <- rep(NA,B)
  
  for(k in 1:B){
    dbk <- datb_cb[[k]]
    tbk <- tb1_cb[[k]]
    theta_hatb[,k] <- (function(theta){f("weibull",dbk,theta,tbk)}) %>% 
      maxLik(.,start=theta_hat,method = "NM",constraints = list(ineqA=matrix(c(1,0,0,1),2,2), ineqB=c(-1e-6,-1e-6))) %>%
      (function(x){x$estimate})
    
    S1_theta_timeb1 <- pweibull(q=as.vector(as.matrix(c(dbk[,1][dbk[,1]>0],dbk[,2][dbk[,2]<Inf]))),shape=theta_hatb[1,k],scale=theta_hatb[2,k],lower.tail=FALSE)
    #S1_theta_timeb2 <- pweibull(q=as.vector(as.matrix(datb_cb[[k]][datb_cb[[k]]$delta2][,1:2])),shape=theta_hatb[1,k],scale=theta_hatb[2,k],lower.tail=FALSE)
    Lamda_nb[k] <- (1-min(tbk)) %>% (function(x){(x-1+tbk)/x}) %>% (function(x){sum((x-S1_theta_timeb1)^2)/2})
  }
  return(list(Lamda_nb,apply(theta_hatb,1,sd)))
}

