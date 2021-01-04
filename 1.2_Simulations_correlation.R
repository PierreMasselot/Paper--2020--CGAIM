###############################################################################
#
#                          Simulation study
#                       Correlation variation
#
###############################################################################

library(parallel)
library(MASS)
library(Matrix)
library(sfsmisc)

# Should be installed from github
# install_github("PierreMasselot/cgaim")
library(cgaim)
source("0_Useful_functions.R")
source("1.0_Benchmark_models.R")

#-------------------------------------------
#     Parameters
#-------------------------------------------

# True index weights alphas
Alpha <- list(
  c(.7, .2, .1, 0), 
  c(0, 0, .5, .5), 
  c(.2, .4, .3, .1)
)

# True ridge functions G
Gfuns <- c(
  function(z, lambda = 1) exp(lambda * scale(z)), 
  function(z, lambda = 5) 1 / (1 + exp(-lambda * scale(z))), 
  function(z, a1 = 50, a2 = .8, b1 = .5, b2 = 2.3) 
    a1 * exp(-b1 * scale(z)) + a2 * exp(b2 * scale(z))
)

# Sample size
n <- 1000

# Number of simulations
ns <- 1000

# Tested correlations
rhovec <- seq(0.25, 0.75, by = 0.25)

#----- Derived objects -----
nr <- length(rhovec)
p <- length(Alpha)
pvec <- sapply(Alpha, length)
ptot <- sum(pvec)

#-------------------------------------------
#     Simulations
#-------------------------------------------

# Initialize cluster for parallel computation
cl <- makeCluster(6)
# Transfer objects in cluster
clusterExport(cl, ls())
clusterEvalQ(cl, {
  library(parallel)
  library(MASS)
  library(Matrix)
  library(cgaim)
})

# Save results
time_samp <- rep(list(vector("numeric", nr)), 4)
g_samp <- rep(list(vector("list", nr)), 4)
z_samp <- rep(list(vector("list", nr)), 4)
alpha_samp <- rep(list(vector("list", nr)), 4)
yhat_samp <- rep(list(vector("list", nr)), 4)

# Loop over simulation designs
for (k in 1:nr){
  print(k); flush.console()
  
  #----- Generate data
  Sigma <- lapply(pvec, function(p){
    sig <- matrix(rhovec[k], nrow = p, ncol = p)
    diag(sig) <- 1
    sig
  })
  X <- Map(MASS::mvrnorm, n = n, mu = lapply(pvec, rep, x = 0), 
    Sigma = Sigma)
  names(X) <- sprintf("X%i", 1:p)
  Xall <- Reduce(cbind, X) # Useful for PPR
  Z <- Map("%*%", X, Alpha)
  G <- mapply(function(g, z) do.call(g, list(z = z)), Gfuns, Z)
  Y <- 5 + rowSums(scale(G))
  Ysim <- replicate(ns, Y + rnorm(n, 0, .2), simplify = F)
  
  # Transfer objects in cluster
  clusterExport(cl, c("Ysim", "X", "Xall", "k", "n"))
  
  #---- CGAIM models
  ## Unconstrained GAIM
  print("GAIM"); flush.console()
  deb <- Sys.time()
  results <- parLapply(cl, Ysim, function(y){
    dat <- c(list(y = y), X)
    cgaim(y ~ g(X1, acons = list(sign.const = 1)) + 
        g(X2, acons = list(sign.const = 1)) + 
        g(X3, acons = list(sign.const = 1)), 
      data = dat, alpha.control = list(norm.type = "sum"))
  })
  time_samp[[1]][k] <- Sys.time() - deb
  # Estimated alphas and functions
  g_samp[[1]][[k]] <- sapply(results, "[[", "gfit", simplify = "array")
  z_samp[[1]][[k]] <- sapply(results, "[[", "indexfit", simplify = "array")
  alphas <- sapply(results, "[[", "alpha")
  alpha_samp[[1]][[k]] <- apply(alphas, 2, unlist)
  yhat_samp[[1]][[k]] <- sapply(results, "[[", "fitted")
  
  ## CGAIM
  print("CGAIM"); flush.console()
  deb <- Sys.time()
  # Apply model
  results <- parLapply(cl, Ysim, function(y){
    dat <- c(list(y = y), X)
    cgaim(y ~ g(X1, fcons = "inc", acons = list(monotone = -1, sign.const = 1)) + 
        g(X2, fcons = "inc", acons = list(monotone = 1, sign.const = 1)) + 
        g(X3, fcons = "cvx", acons = list(sign.const = 1)),
      data = dat, alpha.control = list(norm.type = "sum"),
      smooth.control = list(sp = rep(0, 3)))
  })
  time_samp[[2]][k] <- Sys.time() - deb
  # Estimated alphas and functions
  g_samp[[2]][[k]] <- sapply(results, "[[", "gfit", simplify = "array")
  z_samp[[2]][[k]] <- sapply(results, "[[", "indexfit", simplify = "array")
  alphas <- sapply(results, "[[", "alpha")
  alpha_samp[[2]][[k]] <- apply(alphas, 2, unlist)
  yhat_samp[[2]][[k]] <- sapply(results, "[[", "fitted")
  
  #---- Benchmark models
  ## PPR
  print("PPR"); flush.console()
  deb <- Sys.time()
  results <- parLapply(cl, Ysim, function(y){
    res <- ppr(y = y, x = Xall, nterms = p)  # PPR fitting
    alpha <- Map("[", as.data.frame(res$alpha), split(seq_len(ptot), rep(1:p, pvec)))  # Alphas
    alpha <- Map(cgaim:::normalize, alpha, "sum")
    n <- n
    jf <- 7 + res$smod[1] * (sum(pvec) + 1) # Index for gz
    gz <- matrix(res$smod[jf + 1L:(p * n)], n, p) # gz
    jt <- jf + res$smod[1] * n  # Index for z
    z <- matrix(res$smod[jt + 1L:(p * n)], n, p) # z
    list(alpha = alpha, gz = gz, z = z, coef = c(res$yb, res$beta), 
      fitted = res$fitted.values)
  })
  time_samp[[3]][k] <- Sys.time() - deb
  # Estimated alphas and functions
  g_samp[[3]][[k]] <- sapply(results, "[[", "gz", simplify = "array")
  z_samp[[3]][[k]] <- sapply(results, "[[", "z", simplify = "array")
  alphas <- sapply(results, "[[", "alpha")
  alpha_samp[[3]][[k]] <- apply(alphas, 2, unlist)
  yhat_samp[[3]][[k]] <- sapply(results, "[[", "fitted")
  
  
  ## gMAVE
  print("gMAVE"); flush.console()
  deb <- Sys.time()
  results <- parLapply(cl, Ysim, mave, X, 
    alpha.control = list(norm.type = "sum"))
  time_samp[[4]][k] <- Sys.time() - deb
  # Estimated alphas and functions
  g_samp[[4]][[k]] <- sapply(results, "[[", "gz", simplify = "array")
  z_samp[[4]][[k]] <- sapply(results, "[[", "z", simplify = "array")
  alphas <- sapply(results, "[[", "alpha")
  alpha_samp[[4]][[k]] <- apply(alphas, 2, unlist)
  yhat_samp[[4]][[k]] <- sapply(results, "[[", "fitted")
}

stopCluster(cl)

#---- Save Results
save(rhovec, time_samp, g_samp, alpha_samp, z_samp, yhat_samp, 
  file = "Results/1.2_Simulations_correlation.RData")

