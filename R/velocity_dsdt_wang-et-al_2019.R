vstar <- function(
    fec, # expected per-capita fecundity at N = 0
    disp, # expected dispersal rate at N = 0
    mort, # dispersal mortality
    na.low = FALSE # if TRUE, vstar < 0 are returned as NA, if FALSE, as 0
    ) {
  
  require(assertthat)

  # warning for NA
  if (is.na(fec) | is.na(disp) | is.na(mort)) {
    warning("at least one of the arguments is NA, will return NA")
    return(NA_real_)
  }
  # other tests
  assert_that(
    is.vector(fec), is.vector(disp), is.vector(mort),
    msg = "`fec`, `disp`, and `mort` must all be vectors of length 1"
  )
  assert_that(
    length(fec) == 1, length(disp) == 1, length(mort) == 1,
    msg = "`fec`, `disp`, and `mort` must all be vectors of length 1"
  )
  assert_that(
    is.numeric(fec), is.numeric(disp), is.numeric(mort),
    msg = "`fec`, `disp`, and `mort` must all be numeric"
  )
  assert_that(
    disp >= 0 & disp <= 1,
    msg = "`disp` must be >=0 and <=1"
  )
  assert_that(
    mort >= 0 & mort <= 1,
    msg = "`mort` must be >=0 and <=1"
  )
  assert_that(
    fec >= 0,
    msg = "`fec` must be >=0 (it is per-capita fecundity, not growth rate)"
  )
  assert_that(
    is.logical(na.low), !is.na(na.low),
    msg = "`na.low` must be TRUE or FALSE"
  )


  # formula in Wang et al 2019 (doi: 10.1016/j.tpb.2019.04.003), Appendix A1
  test <- optimise(
    function(kappa,
             rho = fec * (1 - disp * mort), # individuals dying in dispersal are lost to net fecundity
             m = disp * (1 - mort) # and obviously to net dispersal
    ) {
      log(rho * (1 + m * (cosh(kappa) - 1))) / kappa
    },
    lower = 0, upper = 707
    # larger values of upper risk overflow and infinite estimates/NAs
  )
  # Possibly due to the numerical optimisation,
  # test$objective will actually plateau a few 10^-3 or 10^-4 above 1, not 1
  # It will also give negative values if rho < 1
  # so some post-processing needed
  vstar <- test$objective
  if (vstar > 1) {
    vstar <- 1
  }
  if (vstar < 0 & na.low == FALSE) {
    vstar <- 0
  }
  if (vstar < 0 & na.low == TRUE) {
    vstar <- NA_real_
  }
  return(vstar)
}
