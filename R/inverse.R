## The inverse function of the empirical distribution
inverse <- function(u,tb1,x){  #Yt是用来求经验分布的样本 U是随便给的因变量值t
  inn <- !duplicated(x = tb1, fromLast = TRUE)
  tb2 <- tb1[inn]
  x1 <- x[inn]
  tb <- c(1-tb2,1+1e-10)/(1-tb2[length(tb2)])
  V <- rep(NA, length(u))
  tbm <- matrix(c(tb[-(length(tb))],tb[-1]),length(tb)-1,2)
  for(i in 1 : length(u)){
    V[i] <- x1[tbm[,1] <= u[i] & u[i] < tbm[,2]]
  }
  return(V)
}