T1_func <- function(tdist){
  if(tdist=="uniform"){
    T1 <- runif(n = n3, min = min_unif, max = max_unif)
  }else if(tdist=="weibull"){
    T1 <- rweibull(n = n3, shape = shape_weib, scale = scale_weib)
  }else if(tdist=="lognormal"){
    T1 <- rlnorm(n = n3, meanlog  = meanlog, sdlog = sdlog)
  }else if(tdist=="gompertz"){
    T1 <- rgpz(n = n3, lambda = labda, gamma = gama)#1/gama*log(-gama/labda*log(1-runif(n = n3))+1)
  }else if(tdist=="we+5go"){
    u <- runif(n3)
    T1 <- sapply(u, function(x)
      uniroot(function(t, u) pweibull(t,shape_weib,scale_weib,lower.tail=F) - sqrt(5)/10*(pweibull(t,shape_weib,scale_weib,lower.tail=F) - pgpz(t,labda,gama,lower.tail=F))-u, u=x, interval= c(0.1, 100),
              extendInt="yes")$root)
  }else if(tdist=="we+10go"){
    u <- runif(n3)
    T1 <- sapply(u, function(x)
      uniroot(function(t, u) pweibull(t,shape_weib,scale_weib,lower.tail=F) - 1/sqrt(5)*(pweibull(t,shape_weib,scale_weib,lower.tail=F) - pgpz(t,labda,gama,lower.tail=F))-u, u=x, interval= c(0.1, 100),
              extendInt="yes")$root)
  }else if(tdist=="we+15go"){
    u <- runif(n3)
    T1 <- sapply(u, function(x)
      uniroot(function(t, u) pweibull(t,shape_weib,scale_weib,lower.tail=F) - (3*5^(1/2)*(pweibull(t,shape_weib,scale_weib,lower.tail=F) - pgpz(t,labda,gama,lower.tail=F)))/10-u, u=x, interval= c(0.1, 100),
              extendInt="yes")$root)
  }else if(tdist=="we+20go"){
    u <- runif(n3)
    T1 <- sapply(u, function(x)
      uniroot(function(t, u) pweibull(t,shape_weib,scale_weib,lower.tail=F) - (4*5^(1/2)*(pweibull(t,shape_weib,scale_weib,lower.tail=F) - pgpz(t,labda,gama,lower.tail=F)))/10-u, u=x, interval= c(0.1, 100),
              extendInt="yes")$root)
  }else if(tdist=="we+5lo"){
    u <- runif(n3)
    T1 <- sapply(u, function(x)
      uniroot(function(t, u) pweibull(t,shape_weib,scale_weib,lower.tail=F) - (5^(1/2)*(pweibull(t,shape_weib,scale_weib,lower.tail=F) - plnorm(t,meanlog,sdlog,lower.tail=FALSE)))/10-u, u=x, interval= c(0.1, 100),
              extendInt="yes")$root)
  }else if(tdist=="we+10lo"){
    u <- runif(n3)
    T1 <- sapply(u, function(x)
      uniroot(function(t, u) pweibull(t,shape_weib,scale_weib,lower.tail=F) - (5^(1/2)*(pweibull(t,shape_weib,scale_weib,lower.tail=F) - plnorm(t,meanlog,sdlog,lower.tail=FALSE)))/5-u, u=x, interval= c(0.1, 100),
              extendInt="yes")$root)
  }else if(tdist=="we+15lo"){
    u <- runif(n3)
    T1 <- sapply(u, function(x)
      uniroot(function(t, u) pweibull(t,shape_weib,scale_weib,lower.tail=F) - (3*5^(1/2)*(pweibull(t,shape_weib,scale_weib,lower.tail=F) - plnorm(t,meanlog,sdlog,lower.tail=FALSE)))/10-u, u=x, interval= c(0.1, 100),
                extendInt="yes")$root)
  }else if(tdist=="we+20lo"){
    u <- runif(n3)
    T1 <- sapply(u, function(x)
      uniroot(function(t, u) pweibull(t,shape_weib,scale_weib,lower.tail=F) - (4*5^(1/2)*(pweibull(t,shape_weib,scale_weib,lower.tail=F) - plnorm(t,meanlog,sdlog,lower.tail=FALSE)))/10-u, u=x, interval= c(0.1, 100),
              extendInt="yes")$root)
  }
  return(T1)}