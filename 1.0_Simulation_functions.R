#---------------------------
#     Ridge functions
#---------------------------

g1 <- function(z, lambda = 0.5){ 
  exp(lambda * scale(z))
}

g2 <- function(z, lambda = 1){
  1 / (1 + exp(-lambda * scale(z)))
}

g3 <- function(z, a1 = 1, a2 = 1 , b1 = 1, b2 = 1){ 
  a1 * exp(-b1 * scale(z)) + a2 * exp(b2 * scale(z))
}
  
#---------------------------
#  Simulations function
#---------------------------
generate_data <- function(n, ns, Alpha, Gfuns, Gpars, Beta0, Beta1, Xcorr,
  Ysigma)
{
  # Parameters
  p <- length(Alpha)
  pvec <- sapply(Alpha, length)
  d <- sum(pvec)
  pind <- rep(1:p, pvec)
    
  # Simulate X
  Sigma <- lapply(pvec, diag)
  Sigma[[1]][upper.tri(Sigma[[1]])] <- Sigma[[1]][lower.tri(Sigma[[1]])] <- Xcorr
  mu <- lapply(pvec, rep, x = 0)
  X <- Map(MASS::mvrnorm, n = n, mu = mu, Sigma = Sigma)
  names(X) <- sprintf("X%i", 1:p)
  
  # Indices
  Z <- mapply("%*%", X, Alpha)
  
  # Ridge functions
  Gpars <- Map(c, Gpars, z = apply(Z, 2, list))
  G <- mapply(do.call, Gfuns, Gpars)
  G <- scale(G)
  
  # True Y
  Y <- Beta0 + G %*% Beta1 
  
  # Simulate Y with noise
  Ysim <- replicate(ns, Y + Ysigma * rnorm(n), simplify = F)
  return(list(Y = Ysim, X = X, G = G))
}