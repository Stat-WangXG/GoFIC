f <- function(tdist, data, theta, table){
  l <- data$left
  r <- data$right
  dlta1 <- data$delta1
  dlta2 <- data$delta2
  dlta3 <- data$delta3
  phihat <- 1-min(table)
  if(tdist=="uniform"){
    dF1theta <- function(t, theta){
      return(dunif(t, theta[1], theta[2]))}
    S1theta <- function(t, theta){
      return(punif(t, theta[1], theta[2], lower.tail=FALSE))}}
  if(tdist=="weibull"){
    S1theta <- function(t, theta){
      return(pweibull(t,theta[1],theta[2],lower.tail=FALSE))}
    dF1theta <- function(t, theta){
      return(dweibull(t,theta[1],theta[2]))}}
  if(tdist=="lognormal"){
    S1theta <- function(t, theta){
      return(plnorm(t,theta[1],theta[2],lower.tail=FALSE))}
    dF1theta <- function(t, theta){
      return(dlnorm(t,theta[1],theta[2]))}}
  if(tdist=="gompertz"){
    S1theta <- function(t, theta){
      return(pgpz(t,theta[1],theta[2],lower.tail=FALSE))}
      #return(exp(-theta[1]/theta[2]*(exp(theta[2]*t)-1)))}
    dF1theta <- function(t, theta){
      return(dgpz(t,theta[1],theta[2]))}}
      #return(theta[1]*exp(theta[2]*t)*exp(-theta[1]/theta[2]*(exp(theta[2]*t)-1)))}}
  if(tdist=="we+5go"){
    S1theta <- function(t, theta){
      return(pweibull(t,theta[1],theta[2],lower.tail=F) - sqrt(5)/10*(pweibull(t,theta[1],theta[2],lower.tail=F) - pgpz(t,theta[3],theta[4],lower.tail=F)))}
    dF1theta <- function(t, theta){
      return(dweibull(t,theta[1],theta[2]) - sqrt(5)/10*(dweibull(t,theta[1],theta[2]) - dgpz(t,theta[3],theta[4])))}}
  if(tdist=="we+10go"){
    S1theta <- function(t, theta){
      return(pweibull(t,theta[1],theta[2],lower.tail=F) - 1/sqrt(5)*(pweibull(t,theta[1],theta[2],lower.tail=F) - pgpz(t,theta[3],theta[4],lower.tail=F)))}
    dF1theta <- function(t, theta){
      return(dweibull(t,theta[1],theta[2]) - 1/sqrt(5)*(dweibull(t,theta[1],theta[2]) - dgpz(t,theta[3],theta[4])))}}
  if(tdist=="we+15go"){
    S1theta <- function(t, theta){
      return(pweibull(t,theta[1],theta[2],lower.tail=F) - 3*5^(1/2)/10*(pweibull(t,theta[1],theta[2],lower.tail=F) - pgpz(t,theta[3],theta[4],lower.tail=F)))}
    dF1theta <- function(t, theta){
      return(dweibull(t,theta[1],theta[2]) - 3*5^(1/2)/10*(dweibull(t,theta[1],theta[2]) - dgpz(t,theta[3],theta[4])))}}
  if(tdist=="we+20go"){
    S1theta <- function(t, theta){
      return(pweibull(t,theta[1],theta[2],lower.tail=F) - 4*5^(1/2)/10*(pweibull(t,theta[1],theta[2],lower.tail=F) - pgpz(t,theta[3],theta[4],lower.tail=F)))}
    dF1theta <- function(t, theta){
      return(dweibull(t,theta[1],theta[2]) - 4*5^(1/2)/10*(dweibull(t,theta[1],theta[2]) - dgpz(t,theta[3],theta[4])))}}
  if(tdist=="we+5lo"){
    S1theta <- function(t, theta){
      return(pweibull(t,theta[1],theta[2],lower.tail=F) - sqrt(5)/10*(pweibull(t,theta[1],theta[2],lower.tail=F) - plnorm(t,theta[3],theta[4],lower.tail=F)))}
    dF1theta <- function(t, theta){
      return(dweibull(t,theta[1],theta[2]) - sqrt(5)/10*(dweibull(t,theta[1],theta[2]) - dlnorm(t,theta[3],theta[4])))}}
  if(tdist=="we+10lo"){
    S1theta <- function(t, theta){
      return(pweibull(t,theta[1],theta[2],lower.tail=F) - 1/sqrt(5)*(pweibull(t,theta[1],theta[2],lower.tail=F) - plnorm(t,theta[3],theta[4],lower.tail=F)))}
    dF1theta <- function(t, theta){
      return(dweibull(t,theta[1],theta[2]) - 1/sqrt(5)*(dweibull(t,theta[1],theta[2]) - dlnorm(t,theta[3],theta[4])))}}
  if(tdist=="we+15lo"){
    S1theta <- function(t, theta){
      return(pweibull(t,theta[1],theta[2],lower.tail=F) - 3*5^(1/2)/10*(pweibull(t,theta[1],theta[2],lower.tail=F) - plnorm(t,theta[3],theta[4],lower.tail=F)))}
    dF1theta <- function(t, theta){
      return(dweibull(t,theta[1],theta[2]) - 3*5^(1/2)/10*(dweibull(t,theta[1],theta[2]) - dlnorm(t,theta[3],theta[4])))}}
  if(tdist=="we+20lo"){
    S1theta <- function(t, theta){
      return(pweibull(t,theta[1],theta[2],lower.tail=F) - 4*5^(1/2)/10*(pweibull(t,theta[1],theta[2],lower.tail=F) - plnorm(t,theta[3],theta[4],lower.tail=F)))}
    dF1theta <- function(t, theta){
      return(dweibull(t,theta[1],theta[2]) - 4*5^(1/2)/10*(dweibull(t,theta[1],theta[2]) - dlnorm(t,theta[3],theta[4])))}}
  
  deltaS <- log(phihat*(S1theta(l, theta)-S1theta(r, theta)))
  
  logL <- sum(ifelse(dlta1+dlta2+dlta3==0,log(phihat*(dF1theta(l, theta))),0))+
          sum(ifelse(dlta2==0,0,deltaS)+
              ifelse(dlta3==0,0,log(1-phihat+phihat*S1theta(l, theta)))+
              ifelse(dlta1==0,0,log(phihat*(1-S1theta(r, theta)))))
  return (logL)}