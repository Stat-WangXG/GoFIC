pgpz <- function(q,lambda,gamma,lower.tail=TRUE){
  if(lower.tail){
    return(1-exp(-lambda*(exp(gamma*q)-1)/gamma))
  }
  else{
    return(exp(-lambda*(exp(gamma*q)-1)/gamma))
  }
  
}

dgpz <- function(x,lambda,gamma){
  do <- exp(-lambda*(exp(gamma*x)-1)/gamma)
  return(ifelse(do,lambda*exp(gamma*x)*do,0))}

rgpz <- function(n,lambda,gamma){
  log(-gamma*log(runif(n,1e-320,1))/lambda+1)/gamma
}