#' random generator of observation when an effect is expected
#'
#' @param gumb_loc location parameter of extereme distribution
#' @param gumb_scl scale parameter of extereme distribution
#' @param gumb_shp shape parameter of extereme distribution
#'
#' @return one randome value far from zero
#' @export
random_generator_effect <- function(gumb_loc,
                                       gumb_scl,
                                       gumb_shp) {
  fExtremes::rgev(n = 1 ,
                  mu = gumb_loc,
                  beta = gumb_scl,
                  xi = gumb_shp)[1]
}
