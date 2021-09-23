###############################################################
#
#   Implementation of the Groupwise dimention reduction of
#          Li, L., Li, B., Zhu, L.-X., 2010. 
#            Groupwise Dimension Reduction. 
#   Journal of the American Statistical Association
#
###############################################################

# Only estimate for one dimension

gmave <- function(y, x, h, alpha.control = list(), algo.control = list(),
  trace = F)
{
  x <- lapply(x, as.matrix)
  # Initial object
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
  algo.control <- c(algo.control,defalgo.control[!names(defalgo.control) %in% 
    names(algo.control)])
  defalpha.control <- list(norm.type = "L2") #, positive = T)
  alpha.control <- c(alpha.control, 
    defalpha.control[!names(defalpha.control) %in% names(alpha.control)])
  #alpha.control$positive <- rep_len(alpha.control$positive, d)
  #---- Initializing 
  # We consider that dl = 1 for l in 1:M 
  Bvec <- coef(lm(y ~ Xall))[-1]
  Blist <- split(Bvec, pind)
  B <- as.matrix(bdiag(Blist))
  # Prepare Kernel weights
  Kh <- as.matrix(dnorm(dist(Xall) / h) / h) 
  # Prepare difference matrix
  xdiff <- Xall[rep(1:n, n),] - Xall[rep(1:n, each = n),]
  #---- Local linear smoothing step (Step 1)
  predind <- xdiff %*% B
  llrcoef <- matrix(NA, n, p + 1) # Local linear smoothing coefficients
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
    # for (i in 1:n){
    #   invComp <- matrix(0, d, d) # Inverse design matrix
    #   yxmat <- vector("numeric", sum(pvec))# cov matrix
    #   for (j in 1:n){
    #     diffs <- sapply(x, function(xg) xg[j,] - xg[i,])
    #     Rij <- mapply("*", as.data.frame(diffs), llrcoef[i, -1])
    #     Rij <- c(Rij)
    #     invComp <- invComp + (Rij %*% t(Rij) * Kh[i,j])
    #     yxmat <- yxmat + ((y[j] - llrcoef[i,1]) * Kh[i,j] * Rij)
    #   }
    # }
    # iinvComp <- solve(invComp)
    # Bvec2 <- iinvComp %*% yxmat
    Rtot <- xdiff * llrcoefrep[, rep(1:p + 1, pvec)]
    ya <- rep(y, n) - rep(llrcoef[,1], each = n)
    Bvec <- coef(lm(ya ~ 0 + Rtot, weights = c(Kh)))
    # Attempt at coding the positive estimation of Song & Zhu (2016) -> Tout pourrave
#    if (any(alpha.control$positive)){
#      rhos <- cor(cbind(Y,Xall), method = "spearman")[-1,1]
#      lam <- Bvec - apply(iinvComp, 1, sum)
#      lam <- lam * (alpha.control$positive & (sign(Bvec) != sign(rhos)))
#      Lam <- Diagonal(d, lam * sign(rhos))
#      Bvec <- Bvec + iinvComp %*% Lam
#    }
    Blist <- split(Bvec, pind)
    B <- as.matrix(bdiag(Blist)) 
    #---- Local linear smoothing step (Step 1 again)
    predind <- xdiff %*% B
    llrcoef <- matrix(NA, n, p + 1) # Local linear smoothing coefficients
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
  # Normalize and output
  unlist(Map(cgaim:::normalize, Blist, alpha.control$norm.type))
}

###############################################################
#
#           Implementation of the model of
#             Xia, Y., Tong, H., 2005. 
#      Cumulative effects of air pollution on public health. 
#          Statistics in Medicine 25, 3548-3559.  
#
###############################################################

facts <- function(x, y, k = NULL, h = NULL, gshape = "bs", thetashape = "np", 
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
  # isgcons <- shape != "tp"
  # Alist <- vector("list", p)
  # Alist[isgcons] <- lapply(shape[isgcons], function(sh){
  #   sig <- smooth.construct(s(x, bs = sh, k = n), 
  #     data = data.frame(x = 1:n), knots = NULL)$Sigma 
  #   sig <- cbind(1, rbind(0, sig))
  #   bdiag(sig, diag(n))
  # })
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
      # if (isgcons[j] | istcons[j]){
      sij <- singleindex_fit(y = rj, x = x[[j]], k = k, h = h, 
        gshape = gshape[j], B = Blist[[j]], tol = tol, maxit = maxit, 
        dfgrid = dfgrid)
      thetalist[[j]] <- sij$theta
      newgs[,j] <- sij$g
      # } else {
        # bwspec <- npindexbw(ydat = drop(rj), xdat = x[[j]], 
        #   optim.reltol = tol, optim.maxit = maxit)
        # sij <- npindex(bwspec)
        # thetalist[[j]] <- cgaim:::normalize(sij$beta, "L1")
        # newgs[,j] <- sij$mean
      # }
      # Mean centering
      newgs[,j] <- scale(newgs[,j], center = TRUE, scale = FALSE)
    }
    deltaf <- mean((gs - newgs)^2)
    rel.delta <- sqrt(deltaf/fnorm) 
    cbf <- cbf + 1
    gs <- newgs
    fnorm <- mean(gs^2)
    # print(sum((y - beta0 - rowSums(gs))^2))
    # print(rel.delta)
  }
  z <- mapply("%*%", x, thetalist)
  out <- list(theta = thetalist, g = gs, z = z, niter = cbf, intercept = beta0)
}

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
  # ord <- order(z)
  # xord <- x[ord,]
  # yord <- y[ord]
  # Initialize kernel distances
  # xdiff <- xord[rep(1:n, n),] - xord[rep(1:n, each = n),]
  xdiff <- x[rep(1:n, n),] - x[rep(1:n, each = n),]
  xtheta <- xdiff %*% theta
  # kh <- drop(k(xtheta, h(z)))
  kh <- unlist(tapply(k(xtheta, h(z)), rep(1:n, n), function(x) x / sum(x)))
  #----- Step A: Smoothing
  gest <- lapply(dfgrid, function(df) ab_update(y, z, xtheta, gshape, kh, df)) 
  gfun <- gest[[which.min(sapply(gest, "[[", "gcv"))]]$g
  # Initialize MAVE criterion to minimize
  oldobj <- mavecrit(y, xdiff, theta, gfun, kh)
  # Loop
  c1 <- 0
  eps <- tol + 1
  while(eps > tol && c1 <= maxit){
    #----- Step B: update theta
    # Update theta
    theta <- theta_update(y, xdiff, theta, gfun, B, kh)
    # Update index and ordering
    z <- x %*% theta
    # ord <- order(z)
    # xord <- x[ord,]
    # yord <- y[ord]
    # xdiff <- xord[rep(1:n, n),] - xord[rep(1:n, each = n),]
    xdiff <- x[rep(1:n, n),] - x[rep(1:n, each = n),]
    xtheta <- xdiff %*% theta
    kh <- drop(k(xtheta, h(z)))
    #----- Step A: Smoothing
    gest <- lapply(dfgrid, function(df) ab_update(y, z, xtheta, gshape, kh, df)) 
    gfun <- gest[[which.min(sapply(gest, "[[", "gcv"))]]$g
    # Check convergence
    newobj <- mavecrit(y, xdiff, theta, gfun, kh)
    eps <- (oldobj - newobj) / oldobj
    oldobj <- newobj
    c1 <- c1 + 1
    # print(theta)
    # print(eps)
    # print(c1)
  }
  #----- Output
  list(theta = theta, g = gfun[,1], niter = c1)
}

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
  # xthetaDiag <- bdiag(split(xtheta, rep(1:n, each = n)))
  # imat <- bdiag(split(rep(1, n^2), rep(1:n, each = n)))
  # # Design matrix
  # desMat <- cbind(imat, xthetaDiag)
  # Crossproduct matrices for least-squares
  cvec <- crossprod(desMat, kh * rep(y, n))
  dmat <- crossprod(desMat * sqrt(kh))
  if (gshape != "nc"){
    # cvec <- t(A) %*% cvec
    # dmat <- t(A) %*% dmat %*% A
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


# # Previous ab_update code 
# Kh <- dnorm(as.matrix(dist(xord %*% theta, diag = T)) / h) / h
# Dmat <- matrix(0, 2*n, 2*n)
# dvec <- vector("numeric", 2 * n)
# for (i in 1:n){
#   Xmat <- matrix(0, n, 2*n)
#   Xmat[,i] <- 1
#   Xmat[,n + i] <- t(apply(xord, 1, "-", xord[i,])) %*% theta
#   reparaX <- Xmat %*% A
#   wis <- Diagonal(x = Kh[i,])
#   cpres <- crossprod(reparaX, wis %*% reparaX)
#   Dmat <- Dmat + cpres
#   cpres <- crossprod(reparaX, wis %*% yord)
#   dvec <- dvec + cpres
# }
# cMat <- diag(2 * n)[,-1]

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
    # fit <- solve_osqp(2 * Qthe, 2 * Pthe, Cmat, l = rep(0.001, nrow(Cmat)),
    #   pars = osqpSettings(verbose = FALSE, adaptive_rho = FALSE))
    # theta <- B %*% fit$x
    fit <- solve.QP(2 * Qthe, 2 * Pthe, t(Cmat), rep(0.001, nrow(Cmat)))
    theta <- B %*% fit$solution
  }
  # Normalize and return
  cgaim:::normalize(theta, "1")
}

# # Previous theta_update code 
# Qthe <- matrix(0, d, d)
# Pthe <- vector("numeric", d)
# for (i in 1:n){
#   iQ <- matrix(0, d, d)
#   iP <- vector("numeric", d)
#   for (k in 1:n){ 
#     ind <- (i - 1) * n + k
#     dX <- xdiff[ind,]
#     iQ <- iQ + kh[ind] * crossprod(t(dX))
#     iP <- iP + kh[ind] * dX * (y[k] - llrcoefs[i,1])
#   }
#   Qthe <- Qthe + llrcoefs[i,2]^2 * iQ
#   Pthe <- Pthe + llrcoefs[i,2] * iP
# }
# Qthe2 <- t(B) %*% Qthe %*% B
# Pthe2 <- t(B) %*% Pthe

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


###############################################################
#
#           Wrapper for Projection Pursuit Regression
#
###############################################################

# x: list of group matrices
# y: response vector
# norm.type: norm
# ...: parameters to be passed to ppr function
ppr_apply <- function(x, y, norm.type = "L2",  ...)
{
  # Prepare data
  x <- lapply(x, as.matrix)
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
  # Apply ppr
  res <- ppr(y = y, x = Xall, nterms = p, ...)
  # Extract alphas
  alpha <- Map("[", as.data.frame(res$alpha), 
    split(seq_len(ptot), rep(1:p, pvec)))  # Alphas
  unlist(Map(cgaim:::normalize, alpha, norm.type))
  # Extract other elements
  # jf <- 7 + res$smod[1] * (sum(pvec) + 1) # Index for gz
  # gz <- matrix(res$smod[jf + 1L:(p * n)], n, p) # gz
  # jt <- jf + res$smod[1] * n  # Index for z
  # z <- matrix(res$smod[jt + 1L:(p * n)], n, p) # z
}

###############################################################
#
#           Wrapper for SCAIR
#
###############################################################

scair_apply <- function(y, x, shape, norm.type = "L2", ...){
  # Prepare data
  x <- lapply(x, as.matrix)
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
  # Check shape
  if(missing(shape)){
    shape <- rep("in", p)
  }
  # Apply scair
  res <- scair(Xall, y, shape = shape, ...)
  # Extract alphas
  alpha <- Map("[", as.data.frame(res$index), 
    split(seq_len(ptot), rep(1:p, pvec)))  # Alphas
  alpha <- unlist(Map(cgaim:::normalize, alpha, norm.type))
  # Extract betas
  beta <- apply(res$componentfit, 2, sd)
  # Return
  list(alpha = alpha, beta = beta)
}

