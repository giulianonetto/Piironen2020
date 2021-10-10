# Define auxiliary functions ---------------------------------------------------
generate_data_01 <- function(N = 50, p = 500, rho = 0.3, p_rel = p/2) {
  f <- rnorm(n = N, mean = 0, sd = 1)
  y <- rnorm(n = N, mean = f, sd = 1)
  X <- matrix(nrow = N, ncol = p)
  for (i in 1:p) {
    if (i > p_rel) {
      X[, i] <- rnorm(n = N, mean = 0, sd = 1)
    } else {
      X[, i] <- rnorm(n = N, mean = sqrt(rho)*f, sd = sqrt(1 - rho))
    }
  }
  
  return(data.frame(y, X))
}

rmse <- function(pred, truth) {
  sqrt(mean((truth - pred)^2))
}
