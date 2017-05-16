ProblemSolver <- module({
  import("magrittr")
  
  use(.GlobalEnv, attach = TRUE)
  source("ProblemState.R")
  
  start <- function (n, regularizationParam) {
    print("solving started")
    
    #ProblemState$test()
    
    A <- fillMatrix(n)
    b <- fillRightPart(n)
    
    solveSystemWithTihanovReg(matrix=A, rightPart=b, regularizationParam = regularizationParam)
    
    # TODO build solution on inner curve
    # TODO build normal derivetive of solution on inner curve
  }
  
  solveSystemWithTihanovReg <- function(matrix, rightPart, regularizationParam) {
    # regularization parameter
    diagonaMatrix <- diag(dim(matrix)[1])
    matrix <- diagonaMatrix*regularizationParam + t(matrix) %*% matrix
    solve(matrix, rightPart)
  }
  
  initPointsVector <- function (n) {
    step <- pi/n
    seq(0, 2*pi - step, step)
  }
  
  fillMatrix <- function (n) {
    A <- matrix(0, nrow=4*n, ncol=4*n, byrow=TRUE)
    t <- initPointsVector(n)
    for (i in 1:(4*n))
    {
      for (j in 1:(4*n))
      {
        if(i <= 2*n & j <= 2*n){
          A[i,j] <- ProblemState$H11(t[i],t[j]) * pi / n
          next
        }
        
        if(i <= 4*n & j <= 2*n){
          A[i,j] <- ProblemState$H21(t[i- 2 * n],t[j]) * pi / n
          next
        }
        
        if(i <= 2*n & j <= 4*n) {
          A[i,j] <- ProblemState$H12(t[i],t[j - 2 * n]) * pi / n
          if(i == (j - 2 * n)) {
            A[i,j] <- A[i,j] + 0.5/Helpers$vectorNorm(ProblemState$innerCurveDeriv(t[i]))
          }
          next
        }
        
        if(i <= 4*n & j <= 4*n){
          ti <- t[i- 2 * n]
          tj <- t[j - 2 * n]
          A[i,j] <- 2 * pi * ProblemState$H22(ti, tj) * Helpers$hiperSingualtCoeff(ti, n, tj)
          next
        }
      }
    }
    A
  }
  
  fillRightPart <- function (n) {
    t <- initPointsVector(n)
    vec <- rep(0, 4 * n)
    for (i in 1:(2*n)) {
      vec[i] <- ProblemState$solutionOnOuter(t[i])
    }
    
    for (i in (2*n + 1):(4*n)) {
      vec[i] <- ProblemState$solutionDerivetiveOnOuter(t[i - 2 * n])
    }
    vec
  }
})