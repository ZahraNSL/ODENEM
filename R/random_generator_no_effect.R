#' random generator of observation when no effect expected
#'
#' @param trunc_mean mean of truncated normal
#' @param trunc_sd standard deviation of truncated normal
#'
#' @return
#' @export
random_generator_no_effect <- function(trunc_mean,
                                       trunc_sd) {
  abs(msm::rtnorm(1, mean = trunc_mean,
                  sd = trunc_sd))
}
