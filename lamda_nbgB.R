lamda_nbgB <- function(theta_hat,datb_cb,tb1_cb,B){
  theta_hatb <- matrix(NA,2,B)
  Lamda_nb <- rep(NA,B)
  
  for(k in 1:B){
    dbk <- datb_cb[[k]]
    tbk <- tb1_cb[[k]]
    tbk[which(tbk<=0)] <- 1e-30
    theta_hatb[,k] <- (function(theta){f("gompertz",dbk,theta,tbk)}) %>% 
      maxLik(.,start=c(0.1,0.1),method = "NM",constraints = list(ineqA=matrix(c(1,0,0,1),2,2), ineqB=c(-1e-6,-1e-6))) %>%
      (function(x){x$estimate}) #example:start=c(0.1,0.1);simulation:start=theta_hat.
    
    S1_theta_timeb <- pgpz(as.vector(as.matrix(c(dbk[,1][dbk[,1]>0],dbk[,2][dbk[,2]<Inf]))),theta_hatb[1,k],theta_hatb[2,k],lower.tail=F)#exp(-theta_hatb[1,k]/theta_hatb[2,k]*(exp(theta_hatb[2,k]*as.vector(as.matrix(c(dbk[,1][dbk[,1]>0],dbk[,2][dbk[,2]<Inf]))))-1))
    Lamda_nb[k] <- (1-min(tbk)) %>% (function(x){(x-1+tbk)/x}) %>% (function(x){sum((x-S1_theta_timeb)^2)/2})
  }
  return(list(Lamda_nb,apply(theta_hatb,1,sd)))
}

