#'  nonlinear ODE of propagation
#'
#' @param x state
#'
#' @return
#' @export
phi = function(x){
  #1/(1 + exp(-x))
  tanh(x)
}
