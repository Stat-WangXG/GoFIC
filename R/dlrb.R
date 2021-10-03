#==========================delta,left,right=================================
dlrb <- function(Tb1,Tb2,Ub){
  delta1b <- delta2b <- delta3b <- matrix(FALSE, n, B)
  delta1b[(n1+1):n,] <- Tb2 < Ub[1,,]
  delta2b[(n1+1):n,] <- Tb2>Ub[1,,] & Tb2<=Ub[m,,]
  delta3b[(n1+1):n,] <- Tb2 > Ub[m,,]
  
  leftb <- matrix(NA,n,B)
  rightb <- matrix(NA,n,B)
  leftb[delta1b] <- 0
  leftb[delta3b] <- Ub[m,,][which(delta3b[(n1+1):n,])]
  rightb[delta1b] <- Ub[1,,][which(delta1b[(n1+1):n,])]
  rightb[delta3b] <- Inf
  
  for(k in 1:B){
    wh <- which(delta2b[(n1+1):n,k])
    re1 <- Ub[-m,wh,k]
    re2 <- Ub[-1,wh,k]
    re3 <- matrix(rep(Tb2[wh,k],each=m-1),m-1,length(wh))
    wh2 <- which(re3>re1&re3<=re2)
    leftb[wh+n1,k] <- re1[wh2]
    rightb[wh+n1,k] <- re2[wh2]
  }
  leftb[1:n1,] <- rightb[1:n1,] <- Tb1
  
  datb_cb <- split(data.frame(leftb=as.vector(leftb),rightb=as.vector(rightb),delta1b=as.vector(delta1b),delta2b=as.vector(delta2b),delta3b=as.vector(delta3b)),classb)
  
  return(datb_cb)
}
