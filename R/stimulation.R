#' Stimulation function
#'
#' @param stim_val constant value for stimulation function
#' @param t specific time point
#' @param len number of hidden nodes
#'
#' @return decreasing or uprising vector during time
#' @export
stimulation <- function(stim_val,
                        t,
                        len) {
  #return(rep(1 - 1 /(1 + t),len)) # active this line in a case of basic version
  return(rep(stim_val, len))
}
