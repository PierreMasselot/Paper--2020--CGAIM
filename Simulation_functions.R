#---------------------------
#     Ridge functions
#---------------------------

g.lin <- function(z, a = 0, b = 1) a + b * z 
g.cos <- function(z, f = 1) cos(2 * pi * f * z / diff(range(z)))
g.v <- function(z) 1 - exp(-scale(z)^2)
g.sigmoid <- function(z, lambda = 1) 1 / (1 + exp(-lambda * scale(z)))
g.exp <- function(z, lambda = 0.5) exp(lambda * scale(z))
g.jshape <- function(z, a1 = 1, a2 = 1 , b1 = 1, b2 = 1) 
  a1 * exp(-b1 * scale(z)) + a2 * exp(b2 * scale(z)) #  a1 = 50, b1 = .5, b2 = 3
  
#---------------------------
#  Simulations function
#---------------------------
generate_data <- function(n, ns, Alpha, Gfuns, Gpars, Beta0, Beta1, Xcorr,
  Ysigma)
{
  # Derivative parameters and verifications
  p <- length(Alpha)
  pvec <- sapply(Alpha, length)
  d <- sum(pvec)
  pind <- rep(1:p, pvec)
    
  # Simulate X
  Sigma <- lapply(pvec, diag)
  Sigma[[1]][upper.tri(Sigma[[1]])] <- Sigma[[1]][lower.tri(Sigma[[1]])] <- Xcorr
  mu <- lapply(pvec, rep, x = 0)
  X <- Map(mvrnorm, n = n, mu = mu, Sigma = Sigma)
  
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