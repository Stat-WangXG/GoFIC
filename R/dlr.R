dlr <- function(T1,T2,U){
  delta1 <- delta2 <- delta3 <- rep(FALSE,n)
  delta1[(n1+1):n] <- T2 <= U[1,]
  delta2[(n1+1):n] <- (T2>U[1,])&(T2<=U[m,])
  ic <- sum(delta2)/n  #interval censoring
  delta3[(n1+1):n] <- T2 > U[m,]
  rc <- sum(delta3)/n  #right censoring
  
  #============Nonparametric logLikelihood(get phi hat)==============
  left <- rep(NA,n)
  right <- rep(NA,n)
  cS <- colSums(matrix(rep(T2[delta2[(n1+1):n]],each=m), nrow = m, ncol = sum(delta2)) > U[, delta2[(n1+1):n]])
  
  left[delta1] <- 0
  left[delta2] <- U[matrix(c(cS, which(delta2[(n1+1):n])),sum(delta2),2)]
  left[delta3] <- U[m,which(delta3[(n1+1):n])]
  
  right[delta1] <- U[1,which(delta1[(n1+1):n])]
  right[delta2] <- U[matrix(c(cS+1, which(delta2[(n1+1):n])),sum(delta2),2)]
  right[delta3] <- Inf
  
  left[1:n1] <- right[1:n1] <- T1
  dat <- data.frame(left, right, delta1, delta2, delta3)
  return(list(dat,ic,rc))
}