#' bernouli distribution over edge existance probability
#'
#' @param rho Prior over edge existance
#'
#' @return
#' @export
alpha_val = function(rho) {
  return(matrix(
    LaplacesDemon::rbern(length(rho), rho),
    nrow = nrow(rho),
    ncol = ncol(rho)
    #,dimnames =  list(netPara[["S_names"]],netPara[["S_names"]])
  ))
}
