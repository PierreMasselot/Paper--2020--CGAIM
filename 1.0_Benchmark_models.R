###############################################################
#
#   Implementation of the Groupwise dimention reduction of
#          Li, L., Li, B., Zhu, L.-X., 2010. 
#            Groupwise Dimension Reduction. 
# Journal of the American Statistical Association 105, 1188–1201. 
#
###############################################################

mave <- function(y, x, h, alpha.control = list(), algo.control = list()){
  # Initial object
  p <- length(x)
  pvec <- sapply(x, ncol)
  d <- sum(pvec)
  pind <- rep(1:p, pvec)
  Xall <- Reduce(cbind, x)
  n <- nrow(Xall)
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
  B <- bdiag(Blist)
  # Prepare Kernel weights
  Kh <- as.matrix(dnorm(dist(Xall %*% Bvec) / h) / h) 
  #---- Local linear smoothing step (Step 1)
  llrcoef <- matrix(NA, n, p + 1) # Local linear smoothing coefficients
  yhat <- vector("numeric", n)
  for (i in 1:n){    
    dVi <- apply(Xall, 1, "-", Xall[i,])
    Si <- t(B) %*% dVi
    Si <- as.matrix(t(Si))
    llmod <- lm(y ~ Si, weights = Kh[i,])
    llrcoef[i,] <- coef(llmod)
    yhat[i] <- fitted(llmod)[i]
  }
  l2 <- cgaim:::L2(y, yhat)
  eps <- (var(y) - l2)/var(y)
  c1 <- 1
  while(eps > algo.control$tol && c1 <= algo.control$max.iter){    
    #---- Groupwise coefficients updating step (Step 2)
    for (i in 1:n){
      invComp <- matrix(0, d, d) # Inverse design matrix
      yxmat <- vector("numeric", sum(pvec))# cov matrix
      for (j in 1:n){
        diffs <- sapply(x, function(xg) xg[j,] - xg[i,])
        Rij <- mapply("*", diffs, llrcoef[i, -1])
        Rij <- unlist(Rij)
        invComp <- invComp + (Rij %*% t(Rij) * Kh[i,j])
        yxmat <- yxmat + ((y[j] - llrcoef[i,1]) * Kh[i,j] * Rij)
      }      
    }
    iinvComp <- solve(invComp)
    Bvec <- iinvComp %*% yxmat
    # Attempt at coding the positive estimation of Song & Zhu (2016) -> Tout pourrave
#    if (any(alpha.control$positive)){
#      rhos <- cor(cbind(Y,Xall), method = "spearman")[-1,1]
#      lam <- Bvec - apply(iinvComp, 1, sum)
#      lam <- lam * (alpha.control$positive & (sign(Bvec) != sign(rhos)))
#      Lam <- Diagonal(d, lam * sign(rhos))
#      Bvec <- Bvec + iinvComp %*% Lam
#    }
    Blist <- split(Bvec, pind)
    B <- bdiag(Blist) 
    #---- Local linear smoothing step (Step 1 again)
    for (i in 1:n){    
      dVi <- apply(Xall, 1, "-", Xall[i,])
      Si <- t(B) %*% dVi
      llmod <- lm(y ~ as.matrix(t(Si)), weights = Kh[i,])
      llrcoef[i,] <- coef(llmod)
      yhat[i] <- fitted(llmod)[i]
    }
    l2.old <- l2
    l2 <- cgaim:::L2(y, yhat)
    eps <- (l2.old - l2)/l2.old
    c1 <- c1 + 1
#    print(c1)
#    print(eps)
#    flush.console()
  }
  alpha <- Map(cgaim:::normalize, Blist, alpha.control$norm.type)
  z <- mapply("%*%", x, alpha)
  # Rough estimations of the gz functions
  gz <- apply(z, 2, function(zj){
    lowess(zj, llrcoef[,1])$y[rank(zj)]
  })
  gz <- scale(gz)
  betas <- attr(gz, "scaled:scale")
  beta0 <- sum(attr(gz, "scaled:center"))
  out <- list(alpha = alpha, gz = gz, z = z, coef = c(beta0, betas),
    fitted = yhat)
}
