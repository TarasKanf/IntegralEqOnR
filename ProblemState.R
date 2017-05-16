ProblemState <- module({
  import("magrittr")
  
  use(.GlobalEnv, attach = TRUE)
  source("Helpers.R")
  
  test <- function() {
    print("blbabals")
    print(H12(pi/2, 1))
  }
  
  eps <- 0.0001
  innerRadius <- 1
  outerRadius <- 2
  
  #edge values
  solutionOnOuter <- function(t){
    1
  }
  solutionDerivetiveOnOuter <- function (t) {
    0
  }
  
  # curve functions
  innerCurve <- function (t) {
    c(innerRadius * cos(t), innerRadius * sin(t))
  }
  innerCurveDeriv <- function (t) {
    c(-innerRadius * sin(t), innerRadius * cos(t))
  }
  innerCurveDerivSecond <- function (t) {
    c(-innerRadius * cos(t), -innerRadius * sin(t))
  }
  innerCurveNormVec <- function (t) {
    vecNorm <- Helpers$vectorNorm(innerCurve(t))
    deriv <- innerCurveDeriv(t)
    c(deriv[2]/vecNorm, -deriv[1]/vecNorm)
  }
  outerCurve <- function (t) {
    c(outerRadius*cos(t), outerRadius*sin(t))
  }
  outerCurveDeriv <- function (t) {
    c(-outerRadius*sin(t), outerRadius*cos(t))
  }
  outerCurveNormVec <- function (t) {
    vecNorm <- Helpers$vectorNorm(outerCurve(t))
    deriv <- outerCurveDeriv(t)
    c(deriv[2]/vecNorm, -deriv[1]/vecNorm)
  }
  
  # core functions
  H11 <- function (t, tau) {
    curveDeviation <- innerCurve(t) - outerCurve(tau)
    numerator <- Helpers$scalarProduct(curveDeviation, outerCurveNormVec(tau))
    denominator <- Helpers$vectorNorm(curveDeviation)^2
    numerator/denominator/(2*pi)
  }
  
  H12 <- function (t, tau) {
    if(abs(t-tau) > eps){
      curveDeviation <- innerCurve(t) - innerCurve(tau)
      numerator <- Helpers$scalarProduct(curveDeviation, innerCurveNormVec(tau))
      denominator <- Helpers$vectorNorm(curveDeviation)^2
      numerator/denominator/(2*pi)
    } else {
      numerator <- Helpers$scalarProduct(innerCurveDerivSecond(t), innerCurveNormVec(t))
      denominator <- Helpers$vectorNorm(innerCurveDeriv(t))^2
      numerator/denominator/(2*pi)
    }
  }
  
  H21 <- function (t, tau) {
    curveDeviation <- innerCurve(t) - outerCurve(tau)
    numerator1 <- Helpers$scalarProduct(outerCurveNormVec(tau), innerCurveNormVec(t))
    denominator1 <- Helpers$vectorNorm(curveDeviation)^2
    prod21 <- 2 * Helpers$scalarProduct(curveDeviation, outerCurveNormVec(tau))
    numerator2 <- prod21 * Helpers$scalarProduct(curveDeviation, innerCurveNormVec(t))
    denominator2 <- Helpers$vectorNorm(curveDeviation)^4
    
    (numerator1/denominator1 + numerator2/denominator2)/(2*pi)
  }
  
  H22 <- function(t, tau) {
    if(abs(t-tau) > eps){
      curveDeviation <- innerCurve(t) - innerCurve(tau)
      deviationNorm = Helpers$vectorNorm(curveDeviation)
      prod11 <- Helpers$scalarProduct(curveDeviation, innerCurveDeriv(tau))
      prod12 <- Helpers$scalarProduct(curveDeviation, innerCurveDeriv(t))
      numerator1 <- 4 * prod11 * prod12 * sin((t-tau)/2)^2
      denominator1 <- deviationNorm^4
      prod2 <- Helpers$scalarProduct(innerCurveDeriv(tau),innerCurveDeriv(t))
      numerator2 <- prod2 * sin((t-tau)/2)^2
      denominator2 <- deviationNorm^2
      
      (-2 + numerator1/denominator1 + numerator2/denominator2)/(2*pi)
    } else {
      -3/2/(2*pi)
    }
  }
})