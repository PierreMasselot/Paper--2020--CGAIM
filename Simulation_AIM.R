###############################################################################
#
#                          Simulation study
#                       ADDITIVE-INDEX MODELS
#                      Alternative estimations
#
###############################################################################
setwd("C:/Users/masselpl/Documents/Recherche/2017-2019 - Post-doc/Programmes/R/2 - Indices PPR/Paper--2020--CGAIM")

library(parallel)
library(MASS)
library(mgcv)
library(scam)
library(RColorBrewer)
library(fda)
library(forecast)
library(tictoc)
library(gratia)
library(quadprog)
library(GA)
library(scar)
library(osqp)
library(np)
library(devtools)

# Should be installed from github
# install_github("PierreMasselot/cgaim")
library(cgaim)
source("Benchmark_models.R")
source("Simulation_functions.R")


#---- Simulation setting
exp_name <- "Collinearity_study" # Experiment name

constant_parameters <- within(list(),{
  ns <- 5 # Number of simulations
  n <- 1000
  
  Beta0 <- 0  # Intercept of the whole model
  Beta1 <- rep(1 / 3, 3) # Index importance
  Alpha <- list(c(.7, .2, .1), c(.2, .3, .5), c(.2, .6, .2)) # Index weights
  
  Gfuns <- c("g.exp", "g.sigmoid", "g.jshape")  # Ridge function
  Gpars <- list(list(lambda = 1), 
    list(lambda = 5), 
    list(a1 = 50, a2 = 0.8, b1 = .5, b2 = 2.3)
  ) # Parameters for Ridge functions

  Xcorr <- 0

  Ysigma <- .2    # Noise amplitude
})

varying_parameters <- within(list(),{
  Xcorr <- seq(0, 0.8, by = 0.2)  # sample size
})

mod_names <- c("GAIM", "CGAIM", "MGAIM", "PPR", "gMAVE")

#---- Derived parameters
p <- length(constant_parameters$Alpha)
pvec <- sapply(constant_parameters$Alpha, length)
d <- sum(pvec)
pind <- rep(1:p, pvec)
nmod <- length(mod_names)

#-------------------------------------------
#        Loop on varying parameters
#-------------------------------------------

nvp <- length(varying_parameters)
vpvec <- sapply(varying_parameters, length)
vpd <- sum(vpvec)
vpind <- rep(1:nvp, vpvec)
vpval <- unlist(lapply(vpvec, seq))

# Initialize cluster for parallel computation
cl <- makeCluster(2)
# Transfer objects in cluster
clusterExport(cl, ls())
clusterEvalQ(cl, {
  library(parallel)
  library(MASS)
  library(mgcv)
  library(scam)
  library(RColorBrewer)
  library(fda)
  library(forecast)
  library(tictoc)
  library(gratia)
  library(quadprog)
  library(GA)
  library(scar)
  library(osqp)
  library(np)
  library(devtools)
  load_all("C:/Users/masselpl/Documents/Recherche/2017-2019 - Post-doc/Programmes/R/2 - Indices PPR/gaim")
})

# Initialize result objects
datSim <- results <- vector("list", vpd)

exectime <- rep(list(vector("numeric", vpd)), nmod)
g_all <- rep(list(vector("list", vpd)), nmod)
z_all <- rep(list(vector("list", vpd)), nmod)
alpha_all <- rep(list(vector("list", vpd)), nmod)
yhat_all <- rep(list(vector("list", vpd)), nmod)

# Loop over simulation designs
for (k in 1:vpd){
  print(k); flush.console()
  
  #----- Generate data
  kpars <- constant_parameters
  varying_par <- names(varying_parameters)[vpind[k]]
  kpars[[varying_par]] <- varying_parameters[[vpind[k]]][vpval[k]]
  
  datSim[[k]] <- do.call(generate_data, kpars)
  Xall <- Reduce(cbind, datSim[[k]]$X)
  
  # Transfer objects in cluster
  clusterExport(cl, c("Xall", "datSim", "k", "kpars"))
  
  #-------------------------------------------
  #          Apply models
  #-------------------------------------------
  
  #---- Proposed model
  # Unconstrained GAIM
  deb <- Sys.time()
  # Apply model
  results <- parLapply(cl, datSim[[k]]$Y, function(y){
    dat <- data.frame(y, datSim[[k]]$X)
    cgaim(y ~ g(X1, X2, X3) + g(X1.1, X2.1, X3.1) + g(X1.2, X2.2, X3.2),
      data = dat, alpha.control = list(norm.type = "sum"))
  })
  # Execution time
  exectime[[1]][k] <- Sys.time() - deb
  # Estimated alphas and functions
  g_all[[1]][[k]] <- sapply(results, "[[", "gfit", simplify = "array")
  z_all[[1]][[k]] <- sapply(results, "[[", "indexfit", simplify = "array")
  alphas <- sapply(results, "[[", "alpha")
  alpha_all[[1]][[k]] <- apply(alphas, 2, unlist)
  yhat_all[[1]][[k]] <- sapply(results, "[[", "fitted")
  
  # cGAIM
  deb <- Sys.time()
  # Apply model
  results <- parLapply(cl, datSim[[k]]$Y, function(y){
    dat <- data.frame(y, datSim[[k]]$X)
    cgaim(y ~ g(X1, X2, X3, bs = "mpi", 
        constraints = list(monotone = -1, sign.const = 1)) + 
      g(X1.1, X2.1, X3.1, bs = "mpi", 
        constraints = list(monotone = 1, sign.const = 1)) + 
      g(X1.2, X2.2, X3.2, bs = "cx", constraints = list(sign.const = 1)),
      data = dat, alpha.control = list(norm.type = "sum"))
  })
  # Execution time
  exectime[[2]][k] <- Sys.time() - deb
  # Estimated alphas and functions
  g_all[[2]][[k]] <- sapply(results, "[[", "gfit", simplify = "array")
  z_all[[2]][[k]] <- sapply(results, "[[", "indexfit", simplify = "array")
  alphas <- sapply(results, "[[", "alpha")
  alpha_all[[2]][[k]] <- apply(alphas, 2, unlist)
  yhat_all[[2]][[k]] <- sapply(results, "[[", "fitted")
  
  # Misspecified cGAIM
  deb <- Sys.time()
  # Apply model
  results <- parLapply(cl, datSim[[k]]$Y, function(y){
    dat <- data.frame(y, datSim[[k]]$X)
    cgaim(y ~ g(X1, X2, X3, bs = "mpi", 
        constraints = list(sign.const = 1)) + 
      g(X1.1, X2.1, X3.1, bs = "mpi", 
        constraints = list(monotone = -1, sign.const = 1)) + 
      g(X1.2, X2.2, X3.2, bs = "cx", 
        constraints = list(monotone = 1, sign.const = 1)),
      data = dat, alpha.control = list(norm.type = "sum"))
  })
  # Execution time
  exectime[[3]][k] <- Sys.time() - deb
  # Estimated alphas and functions
  g_all[[3]][[k]] <- sapply(results, "[[", "gfit", simplify = "array")
  z_all[[3]][[k]] <- sapply(results, "[[", "indexfit", simplify = "array")
  alphas <- sapply(results, "[[", "alpha")
  alpha_all[[3]][[k]] <- apply(alphas, 2, unlist)
  yhat_all[[3]][[k]] <- sapply(results, "[[", "fitted")  
  
  #---- Benchmark models
  # PPR
  deb <- Sys.time()
  # Apply model
  results <- parLapply(cl, datSim[[k]]$Y, function(y){
    res <- ppr(y = y, x = Xall, nterms = p)  # PPR fitting
    alpha <- Map("[", as.data.frame(res$alpha), split(1:9, pind))  # Alphas
    alpha <- Map(normalize, alpha, "sum")
    n <- kpars$n
    jf <- 7 + res$smod[1] * (sum(pvec) + 1) # Index for gz
    gz <- matrix(res$smod[jf + 1L:(p * n)], n, p) # gz
    jt <- jf + res$smod[1] * n  # Index for z
    z <- matrix(res$smod[jt + 1L:(p * n)], n, p) # z
    list(alpha = alpha, gz = gz, z = z, coef = c(res$yb, res$beta), 
      fitted = res$fitted.values)
  })
  # Execution time
  exectime[[4]][k] <- Sys.time() - deb
  # Estimated alphas and functions
  g_all[[4]][[k]] <- sapply(results, "[[", "gz", simplify = "array")
  z_all[[4]][[k]] <- sapply(results, "[[", "z", simplify = "array")
  alphas <- sapply(results, "[[", "alpha")
  alpha_all[[4]][[k]] <- apply(alphas, 2, unlist)
  yhat_all[[4]][[k]] <- sapply(results, "[[", "fitted")
  
  
  # gMAVE (Li et al. 2010)
  deb <- Sys.time()
  # Apply model
  results <- parLapply(cl, datSim[[k]]$Y, mave, datSim[[k]]$X, 
    alpha.control = list(norm.type = "sum"))
  # Execution time
  exectime[[5]][k] <- Sys.time() - deb
  # Estimated alphas and functions
  g_all[[5]][[k]] <- sapply(results, "[[", "gz", simplify = "array")
  z_all[[5]][[k]] <- sapply(results, "[[", "z", simplify = "array")
  alphas <- sapply(results, "[[", "alpha")
  alpha_all[[5]][[k]] <- apply(alphas, 2, unlist)
  yhat_all[[5]][[k]] <- sapply(results, "[[", "fitted")
}

stopCluster(cl)

#---- Save Results
save(constant_parameters, varying_parameters, datSim, exectime, g_all, 
  alpha_all, z_all, yhat_all, file = sprintf("%s_resultsV3.RData", exp_name))

