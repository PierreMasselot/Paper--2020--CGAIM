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
  defalgo.control <- list(tol = 5e-3, max.iter = 50)
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
  l2 <- L2(y, yhat)
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
    l2 <- L2(y, yhat)
    eps <- (l2.old - l2)/l2.old
    c1 <- c1 + 1
#    print(c1)
#    print(eps)
#    flush.console()
  }
  alpha <- Map(normalize, Blist, alpha.control$norm.type)
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

###############################################################
#
#           Implementation of the model of
#             Xia, Y., Tong, H., 2005. 
#      Cumulative effects of air pollution on public health. 
#          Statistics in Medicine 25, 3548–3559.  
#
###############################################################

backfit_llr <- function(y, x, h, alpha.control = list(), 
  algo.control = list(), smooth.control = list())
{
  n <- length(y)
  p <- length(x)
  pvec <- sapply(x, ncol)
  d <- sum(pvec)
  if (missing(h)) h <- ((4/(d+2))^(1/(d+4))) * (n^(-1/(d+4)))
  # Default algorithm control
  def.control <- list(bf.tol = 5e-03, bf.maxit = 50, tol = 5e-3, max.iter = 50)
  algo.control <- c(algo.control, 
    def.control[!names(def.control) %in% names(algo.control)])
  defalpha.control <- list(norm.type = "L2", monotone = 0, sign.const = 0)
  alpha.control <- c(alpha.control, 
    defalpha.control[!names(defalpha.control) %in% names(alpha.control)])
  alpha.control[c("norm.type", "sign.const", "monotone")] <- 
    lapply(alpha.control[c("norm.type", "sign.const", "monotone")], 
    rep_len, length.out = p)
  defsmoo.control <- list(shape = "")
  smooth.control <- c(smooth.control, 
    defsmoo.control[!names(defsmoo.control) %in% names(smooth.control)])
  smooth.control$shape <- rep_len(smooth.control$shape, length.out = p)
  # Initialization
  gs <- matrix(0, n, p)
  rel.delta <- algo.control$bf.tol + 1
  cbf <- 0
  beta0 <- mean(y)
  y <- y - beta0
  alpha <- vector("list", p)
  z <- matrix(NA, n, p)
  # Backfitting
  while (rel.delta > algo.control$bf.tol && cbf < algo.control$bf.maxit){
    deltaf <- 0
    f0norm <- sum(apply(gs^2, 2, mean))
    for (j in 1:p){
      rj <- y - rowSums(gs[,-j, drop = F])
      if (alpha.control$monotone[j] == 0 && smooth.control$shape == ""){
        bwspec <- npindexbw(ydat = as.vector(rj), xdat = x[[j]])
        fit <- npindex(bwspec)
        alpha[[j]] <- fit$beta
        z[,j] <- fit$index
        deltaf <- deltaf + mean((fit$mean - gs[,j])^2)
        gs[,j] <- fit$mean
      } else {
        theta <- (1:pvec[j] - 1) / sum((1:pvec[j] - 1))
        if (smooth.control$shape[j] == ""){
          Sigmat <- diag(n)
        } else {
          Sigmat <- smooth.construct(s(V1, bs = smooth.control$shape[j], k = n), 
            data = as.data.frame(x = rj), knots = NULL)$Sigma
          Sigmat <- cbind(1, rbind(0, Sigmat))
        }
        A <- bdiag(Sigmat, diag(n))
        alphaCons <- switch(alpha.control$monotone[j] + 2, "mpd", "tp", "mpi")
        if (alphaCons == "tp"){
          B <- diag(pvec[j])
        } else {
          B <- smooth.construct(s(V1, bs = alphaCons, k = max(5, pvec[j])), 
          data = as.data.frame(x = rj), knots = NULL)$Sigma
          B <- B[1:(pvec[j]-1), 1:(pvec[j]-1)]
          B <- cbind(1, rbind(0, B))
        }
        #----- Step A (first time)
        zj <- x[[j]] %*% theta
        xord <- x[[j]][order(zj),]
        rjord <- rj[order(zj)]
        yord <- y[order(zj)]
        Kh <- dnorm(as.matrix(dist(xord %*% theta, diag = T)) / h) / h
        Dmat <- matrix(0, 2*n, 2*n)
        dvec <- vector("numeric", 2 * n)
        for (i in 1:n){
          Xmat <- matrix(0, n, 2*n)
          Xmat[,i] <- 1
          Xmat[,n + i] <- t(apply(xord, 1, "-", xord[i,])) %*% theta
          reparaX <- Xmat %*% A
          wis <- Diagonal(x = Kh[i,]) 
          cpres <- crossprod(reparaX, wis %*% reparaX)
          Dmat <- Dmat + cpres
          cpres <- crossprod(reparaX, wis %*% yord)
          dvec <- dvec + cpres        
        }
        cMat <- diag(2 * n)[,-1]        
        fit <- solve.QP(Dmat, dvec, cMat)
        llrcoefs <- A %*% fit$solution
        llrcoefs <- matrix(llrcoefs, ncol = 2)
        l2 <- L2(yord, llrcoefs[,1])     
        c1 <- 0
        eps <- algo.control$tol + 1
        while(eps > algo.control$tol && c1 <= algo.control$max.iter){           
          #----- Step B
          Qthe <- matrix(0, pvec[j], pvec[j])
          Pthe <- vector("numeric", pvec[j])
          for (i in 1:n){
            iQ <- matrix(0, pvec[j], pvec[j])
            iP <- vector("numeric", pvec[j])
            for (k in 1:n){            
              dX <- xord[k,] - xord[i,]
              iQ <- iQ + Kh[k,i] * crossprod(t(dX))
              iP <- iP + Kh[k,i] * dX * (rjord[k] - llrcoefs[i,1])
            }
            Qthe <- Qthe + llrcoefs[i,2]^2 * iQ
            Pthe <- Pthe + llrcoefs[i,2] * iP
          }
          Qthe <- t(B) %*% Qthe %*% B
          Pthe <- t(B) %*% Pthe
          fit <- solve.QP(Qthe, Pthe, diag(pvec[j])) 
          theta <- B %*% fit$solution
          theta <- normalize(theta, alpha.control$norm.type[j])
          # ---- Step A  (again)
          zj <- x[[j]] %*% theta
          xord <- x[[j]][order(zj),]
          rjord <- rj[order(zj)]
          yord <- y[order(zj)]
          Kh <- dnorm(as.matrix(dist(xord %*% theta, diag = T)) / h) / h
          Dmat <- matrix(0, 2*n, 2*n)
          dvec <- vector("numeric", 2 * n)
          for (i in 1:n){
            Xmat <- matrix(0, n, 2*n)
            Xmat[,i] <- 1
            Xmat[,n + i] <- t(apply(xord, 1, "-", xord[i,])) %*% theta
            reparaX <- Xmat %*% A
            wis <- Diagonal(x = Kh[i,]) 
            cpres <- crossprod(reparaX, wis %*% reparaX)
            Dmat <- Dmat + cpres
            cpres <- crossprod(reparaX, wis %*% yord)
            dvec <- dvec + cpres        
          }
          cMat <- diag(2 * n)[,-1]        
          fit <- solve.QP(Dmat, dvec, cMat)
          llrcoefs <- A %*% fit$solution
          llrcoefs <- matrix(llrcoefs, ncol = 2)
          l2.old <- l2     
          l2 <- L2(yord, llrcoefs[,1])
          eps <- (l2.old - l2)/l2.old
          c1 <- c1 + 1 
        }
        alpha[[j]] <- theta
        z[,j] <- zj
        deltaf <- deltaf + mean((llrcoefs[rank(zj),1] - gs[,j])^2)
        gs[,j] <- llrcoefs[rank(zj),1] 
        alpha[[j]] <- theta
      } 
    }
    rel.delta <- sqrt(deltaf/f0norm) 
    cbf <- cbf + 1
  }
  gz <- scale(gs)  # gz is a rough approximation
  betas <- attr(gz, "scaled:scale")
  beta0 <- beta0 + sum(attr(gz, "scaled:center"))
  out <- list(alpha = alpha, gz = gz, z = z, coef = c(beta0, betas),
    fitted = beta0 + gz %*% betas)
}