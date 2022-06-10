################################################################################
#
#   R code implementing the FACTS model used in simulation in
#
#   Masselot et al., 2022
#   Constrained groupwise additive index models
#   Biostatistics
#
#   Original method:
#   Kong et al., 2010. 
#   Statistical modelling of nonlinear long-term cumulative effects. 
#   Statistica Sinica
#
#   Author: Pierre Masselot
#
################################################################################

#------------------------
# Main function
#------------------------

# The algorithm proposed by Kong et al. (2010) is a backfitting which loops
# across indices until the fit does not progresses. It iteratively fits single 
# index models through local linear estimation.

#----- Parameters
# x: A list of matrices representing the groups of predictors.
# y: The response vector
# k: A kernel function for smoothing. Default to the Gaussian kernel.
# h: Bandwidth for kernel smoothing. Default to the Silverman's rule of thumb.
# gshape: The shape constraint to impose for each index. 
#     One of 'nc' (no constraint, the default), 'inc' increasing,
#     'dec' (decreasing), 'cvx' (convex) or 'ccv' (concave).
# thetashape: constraint on group weights theta. Either 'nc'
#     (no constraint) or 'inc' (increasing).
# dfgrid: grid of degrees of freedom for the constrained splines used for
#     shape-constrained smoothing. At each step of the algorithm the best
#     degrees of freedom are chosen by GCV.
# bf.tol: tolerance for convergence of the MAVE criterion in the backfitting
#     algorithm.
# bf.maxit: maximum number of iterations in the backfitting algorithm.
# tol: tolerance for convergence of algorithm in each single index fitting.
# maxit: maximum number of iterations for sinlge index fitting.

#----- Value
# A list with elements
# theta: list of estimated index weights.
# g: estimated ridge functions.
# z: evaluated indices using estimated thetas.
# niter: number of backfitting iteration performed.
# intercept: estimated intercept of the model.

facts <- function(x, y, k = NULL, h = NULL, gshape = "nc", thetashape = "nc", 
  dfgrid = 5:10, bf.tol = 1e-3, bf.maxit = 50, tol = 1e-3, maxit = 50)
{
  # Prepare data
  x <- lapply(x, as.matrix)
  y <- as.vector(y)
  n <- unique(c(length(y), sapply(x, nrow)))
  if (length(n) != 1) stop("inconsistent number of observations in x and y")
  # Dimensions
  p <- length(x)
  pvec <- sapply(x, ncol)
  d <- sum(pvec)
  # Prepare kernel function and bandwidth
  if (is.null(k)) k <- function(x, h) dnorm(x / h) / h
  if (is.null(h)) h <- function(x) 1.06 * sd(x) / (length(x) ^ (1/5))
  # Prepare smoothing constraints
  gshape <- match.arg(gshape, c("nc", "inc", "dec", "cvx", "ccv"), 
    several.ok = T)
  gshape <- rep_len(gshape, p)
  # Prepare index coefficient constraints  
  thetashape <- match.arg(thetashape, c("nc", "inc"), 
    several.ok = T)
  thetashape <- rep_len(thetashape, p)
  istcons <- thetashape == "inc"
  Blist <- vector("list", p)
  Blist[istcons] <- Map(function(sh, d){
    sig <- smooth.construct(s(x, bs = "mpi", k = d), 
      data = data.frame(x = 1:n), knots = NULL)$Sigma 
    cbind(1, rbind(0, sig))
  }, thetashape[istcons], pvec[istcons])  
  # Initialization
  gs <- newgs <- matrix(0, n, p)
  thetalist <- vector("list", p)
  fnorm <- 1
  rel.delta <- bf.tol + 1
  cbf <- 0
  beta0 <- mean(y)
  # Backfitting
  while (rel.delta > bf.tol && cbf < bf.maxit){
    for (j in 1:p){
      rj <- y - beta0 - rowSums(gs[,-j, drop = F])
      # Single index fitting
      sij <- singleindex_fit(y = rj, x = x[[j]], k = k, h = h, 
        gshape = gshape[j], B = Blist[[j]], tol = tol, maxit = maxit, 
        dfgrid = dfgrid)
      thetalist[[j]] <- sij$theta
      newgs[,j] <- sij$g
      # Mean centering
      newgs[,j] <- scale(newgs[,j], center = TRUE, scale = FALSE)
    }
    deltaf <- mean((gs - newgs)^2)
    rel.delta <- sqrt(deltaf/fnorm) 
    cbf <- cbf + 1
    gs <- newgs
    fnorm <- mean(gs^2)
  }
  z <- mapply("%*%", x, thetalist)
  out <- list(theta = thetalist, g = gs, z = z, niter = cbf, intercept = beta0)
}

#------------------------
# Internal functions
#------------------------

#----- Fits single index model within each iteration of the backfitting algorithm

singleindex_fit <- function(y, x, k, h, gshape, B, 
  tol = 1e-3, maxit = 50, dfgrid)
{
  d <- ncol(x)
  n <- nrow(x)
  # Initialize theta
  theta <- rep(1/d, d)
  if (!is.null(B)) theta <- B %*% theta
  # Compute z and order
  z <- x %*% theta
  # Initialize kernel distances
  xdiff <- x[rep(1:n, n),] - x[rep(1:n, each = n),]
  xtheta <- xdiff %*% theta
  kh <- unlist(tapply(k(xtheta, h(z)), rep(1:n, n), function(x) x / sum(x)))
  # Step A: Smoothing
  gest <- lapply(dfgrid, function(df) ab_update(y, z, xtheta, gshape, kh, df)) 
  gfun <- gest[[which.min(sapply(gest, "[[", "gcv"))]]$g
  # Initialize MAVE criterion to minimize
  oldobj <- mavecrit(y, xdiff, theta, gfun, kh)
  # Loop
  c1 <- 0
  eps <- tol + 1
  while(eps > tol && c1 <= maxit){
    # Step B: update theta
    theta <- theta_update(y, xdiff, theta, gfun, B, kh)
    # Update index and ordering
    z <- x %*% theta
    xdiff <- x[rep(1:n, n),] - x[rep(1:n, each = n),]
    xtheta <- xdiff %*% theta
    kh <- drop(k(xtheta, h(z)))
    # Step A: Smoothing
    gest <- lapply(dfgrid, function(df) ab_update(y, z, xtheta, gshape, kh, df)) 
    gfun <- gest[[which.min(sapply(gest, "[[", "gcv"))]]$g
    # Check convergence
    newobj <- mavecrit(y, xdiff, theta, gfun, kh)
    eps <- (oldobj - newobj) / oldobj
    oldobj <- newobj
    c1 <- c1 + 1
  }
  # Output
  list(theta = theta, g = gfun[,1], niter = c1)
}

#----- Update smooth function
ab_update <- function(y, z, xtheta, gshape, kh, df)
{
  # Dimensions
  n <- length(y)
  # Create design matrix of spline smoothing
  splMat <- if(gshape %in% c("inc", "dec")){
    cbind(1, iSpline(z, df = df, intercept = T))
  } else if (gshape %in% c("cvx", "ccv")){
    cbind(1, cSpline(z, df = df, intercept = T))
  } else {
    bSpline(z, df = df, intercept = T)
  }
  dsplMat <- if(gshape %in% c("inc", "dec")){
    cbind(0, iSpline(z, df = df, intercept = T, derivs = 1))
  } else if (gshape %in% c("cvx", "ccv")){
    cbind(0, cSpline(z, df = df, intercept = T, derivs = 1))
  } else {
    bSpline(z, df = df, intercept = T, derivs = 1)
  }
  desMat <- splMat[rep(1:n, n),] + dsplMat[rep(1:n, n),] * drop(xtheta)
  # Crossproduct matrices for least-squares
  cvec <- crossprod(desMat, kh * rep(y, n))
  dmat <- crossprod(desMat * sqrt(kh))
  if (gshape != "nc"){
    Cmat <- cbind(0, diag(df))
    if (gshape %in% c("dec", "ccv")) Cmat <- -Cmat
    # Solve QP
    fit <- solve.QP(2 * dmat, 2 * cvec, t(Cmat), bvec = rep(0.001, nrow(Cmat)))
    gres <- splMat %*% fit$solution
    dgres <- dsplMat %*% fit$solution
  } else {
    fit <- solve(dmat) %*% cvec
    gres <- splMat %*% fit
    dgres <- dsplMat %*% fit
  }
  # Compute GCV
  hatmat <- splMat %*% solve(crossprod(splMat)) %*% t(splMat)
  gcv <- mean((y - gres)^2) / (1 - sum(diag(hatmat)) / n)^2
  # Output
  list(g = cbind(gres, dgres), gcv = gcv)
}

#----- Update index weights theta
theta_update <- function(y, xdiff, theta, llrcoefs, B, kh)
{
  d <- ncol(xdiff)
  n <- length(y)
  # Compute Q matrix
  llrexp <- llrcoefs[rep(1:n, each = n),]
  xouter <- xdiff[,rep(1:d, d)] * xdiff[,rep(1:d, each = d)] * kh * 
    (llrexp[,2] ^ 2)
  Qthe <- matrix(colSums(xouter), d, d)
  # Compute P matrix
  Pthe <- colSums(xdiff * kh * (rep(y, n) - llrexp[,1]) * llrexp[,2])
  # Estimate coefficients
  if (is.null(B)){
    theta <- solve(Qthe) %*% Pthe
  } else {
    Qthe <- t(B) %*% Qthe %*% B
    Pthe <- t(B) %*% Pthe
    Cmat <- diag(d)
    fit <- solve.QP(2 * Qthe, 2 * Pthe, t(Cmat), rep(0.001, nrow(Cmat)))
    theta <- B %*% fit$solution
  }
  # Normalize and return
  cgaim:::normalize(theta, "1")
}

#----- Compute MAVE criterion
mavecrit <- function(y, xdiff, theta, llrcoefs, kh)
{
  n <- length(y)
  # Design matrix
  Xtheta <- xdiff %*% theta
  Xdes <- cbind(1, Xtheta)
  # Compute squared deviations
  fits <- rowSums(Xdes * llrcoefs[rep(1:n, each = n),])
  err <- kh * (y - fits)^2
  # return mave criterion
  sum(err)
}
