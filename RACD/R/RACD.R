#'@import Rcpp
#'@import RcppProgress
#'@importFrom deSolve ode
#'@useDynLib RACD
NULL

fequal <- function(x,y,tol=.Machine$double.eps^2){
  abs(x-y) < tol
}
