getsurv2 <- function(data){
  L <- data$left
  R <- data$right
  x <- c(L[L>0],R[R<Inf])
  g <- getsurv(times = x, icfit=icfit(L,R))[[1]]$S
  return(g)
}