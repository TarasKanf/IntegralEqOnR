Helpers <- module({
  use(.GlobalEnv, attach = TRUE)
  
  scalarProduct <- function (vec1, vec2) {
    sum(vec1*vec2)
  }
  
  vectorNorm <- function (vec) {
    sqrt(sum(vec*vec))
  }
  
  hiperSingularCoeff <- function (t, n, tj) {
    sum <- 0
    delta_t <- t - tj
    for (i in 1:(n-1))
    {
      sum = sum + i * cos(i * delta_t)
    }
    sum <- sum * (-1 / n)
    sum <- sum - cos(n * delta_t) / 2
  }
})
