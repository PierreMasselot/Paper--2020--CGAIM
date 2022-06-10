################################################################################
#
#   R code implementing the gMAVE model used in simulation in
#
#   Masselot et al., 2022
#   Constrained groupwise additive index models
#   Biostatistics
#
#   Original method:
#   Li, et al., 2010. 
#   Groupwise Dimension Reduction. 
#   Journal of the American Statistical Association
#
#   Author: Pierre Masselot
#
################################################################################

#----- Parameters
# y: The response vector
# x: A list of matrices representing the groups of predictors.
# h: The bandwidth for local-linear smoothing. If missing, computed automatically
#     using the formula given in Li et al. (2010) JASA
# alpha.norm: The norm constraint for each alpha vector. See cgaim package for
#     list of possible norms.
# algo.control: A list of optional parameters controlling the algorithm.
#     'tol' the tolerance on the gMAVE criterion for convergence (default to 0.001) 
#     'max.iter' maximum number of iterations (default to 50)
# trace: Logical. If true, prints evolution of number of iteration and gMAVE
#     criterion.

#----- Value
# A vector of stacked weights alpha. Follows the order of variables in x.

gmave <- function(y, x, h, alpha.norm = "L2", algo.control = list(),
  trace = F)
{
  x <- lapply(x, as.matrix)
  # Initialize objects
  p <- length(x)
  pvec <- sapply(x, ncol)
  d <- sum(pvec)
  pind <- rep(1:p, pvec)
  Xall <- Reduce(cbind, x)
  nx <- sapply(x, nrow)
  ny <- length(y)
  n <- unique(c(nx, ny))
  if (length(n) != 1){
    stop("inconsistent number of observation")
  }
  # Default parameters for the algorithm
  if (missing(h)) h <- ((4/(d+2))^(1/(d+4))) * (n^(-1/(d+4)))
  defalgo.control <- list(tol = 1e-3, max.iter = 50)
  algo.control <- c(algo.control, defalgo.control[!names(defalgo.control) %in% 
    names(algo.control)])
  # Initializing alphas
  Bvec <- coef(lm(y ~ Xall))[-1]
  Blist <- split(Bvec, pind)
  B <- as.matrix(bdiag(Blist))
  # Prepare Kernel weights
  Kh <- as.matrix(dnorm(dist(Xall) / h) / h) 
  # Prepare difference matrix
  xdiff <- Xall[rep(1:n, n),] - Xall[rep(1:n, each = n),]
  #---- Local linear smoothing step (Step 1)
  predind <- xdiff %*% B
  llrcoef <- matrix(NA, n, p + 1)
  for (i in 1:n){    
    Si <- predind[(i - 1) * n + 1:n,]
    llmod <- lm(y ~ Si, weights = Kh[i,])
    llrcoef[i,] <- coef(llmod)
  }
  llrcoefrep <- llrcoef[rep(1:n, each = n),]
  allpreds <- rowSums(cbind(1, predind) * llrcoefrep)
  gmavecritold <- sum(c(Kh) * (y[rep(1:n, n)] - allpreds)^2)
  # Initialize convergence criteria
  eps <- 1
  c1 <- 0
  # Trace
  if (trace){
    print(c1)
    print(gmavecrit)
    flush.console()
  }
  # Loop
  while(eps > algo.control$tol && c1 <= algo.control$max.iter){    
    #---- Groupwise coefficients updating step (Step 2)
    Rtot <- xdiff * llrcoefrep[, rep(1:p + 1, pvec)]
    ya <- rep(y, n) - rep(llrcoef[,1], each = n)
    Bvec <- coef(lm(ya ~ 0 + Rtot, weights = c(Kh)))
    Blist <- split(Bvec, pind)
    B <- as.matrix(bdiag(Blist)) 
    #---- Local linear smoothing step (Step 1 again)
    predind <- xdiff %*% B
    llrcoef <- matrix(NA, n, p + 1)
    for (i in 1:n){    
      Si <- predind[(i - 1) * n + 1:n,]
      llmod <- lm(y ~ Si, weights = Kh[i,])
      llrcoef[i,] <- coef(llmod)
    }
    # Compute gmave criterion for convergence
    allpreds <- rowSums(cbind(1, predind) * llrcoef[rep(1:n, each = n),])
    gmavecrit <- sum(c(Kh) * (y[rep(1:n, n)] - allpreds)^2)
    # Update
    eps <- (gmavecritold - gmavecrit) / gmavecritold
    c1 <- c1 + 1
    gmavecritold <- gmavecrit
    # Trace
    if(trace){
      print(c1)
      print(gmavecrit)
      flush.console()
    }
  }
  # Normalize and return
  unlist(Map(cgaim:::normalize, Blist, alpha.norm))
}
