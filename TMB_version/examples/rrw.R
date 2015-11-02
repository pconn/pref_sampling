
#' SIMULATE AN ICAR PROCESS
#' @param Q Precision matrix for the ICAR process
#' @return Spatial random effects
#' @export
#' @keywords ICAR, simulation
#' @author Devin Johnson
rrw <- function(Q){
  v <- eigen(Q, TRUE)
  val.inv <- sqrt(ifelse(v$values>sqrt(.Machine$double.eps), 1/v$values, 0))
  P <- v$vectors
  sim <- P%*%diag(val.inv)%*%rnorm(dim(Q)[1], 0, 1)
  X <- rep(1,length(sim))
  if(sum(val.inv==0)==2) X <- cbind(X, 1:length(sim))
  sim <- sim-X%*%solve(crossprod(X), crossprod(X,sim))
  return(sim)
}
