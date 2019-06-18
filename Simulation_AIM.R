###############################################################################
#
#                          Simulation study
#                       ADDITIVE-INDEX MODELS
#                      Alternative estimations
#
###############################################################################
setwd("D:/Pierre Masselot/AIM")

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

source("gaim/Secondary functions.R")
source("gaim/aim.R")
source("Benchmark_models.R")


#---- Functions
# Ridge function
g.lin <- function(z, a = 0, b = 1) a + b * z 
g.cos <- function(z, f = 1) cos(2 * pi * f * z / diff(range(z)))
g.v <- function(z) 1 - exp(-scale(z)^2)
g.sigmoid <- function(z, lambda = 1) 1 / (1 + exp(-lambda * scale(z)))
g.exp <- function(z, lambda = 0.5) exp(lambda * scale(z))
g.jshape <- function(z, a1 = 1, a2 = 1 , b1 = 1, b2 = 1) 
  a1 * exp(-b1 * scale(z)) + a2 * exp(b2 * scale(z)) #  a1 = 50, b1 = .5, b2 = 3

# Covariance matrix
sig.dist <- function(p, sig = .5) 
  sig ^ as.matrix(dist(1:p, upper = T, diag = T)) 
sig.identity <- function(p) diag(p)

#---- Simulation setting
parameters <- within(list(),{
  ns <- 1000 # Number of simulations
  n <- 1000 # sample size
  
  Beta0 <- 0  # Intercept of the whole model
  Beta1 <- rep(1 / 3, 3) # Index importance
  Alpha <- list(c(.7, .2, .1), c(.1, .9), c(.1, .4, .3, .2)) # Index weights
  
  Gfuns <- c("g.exp", "g.sigmoid", "g.jshape")  # Ridge function
  Gpars <- list(list(lambda = 1), 
    list(lambda = 5), 
    list(a1 = 50, a2 = 0.8, b1 = .5, b2 = 2.3)
  ) # Parameters for Ridge functions
  
  XSfun <- c("sig.dist", "sig.dist", "sig.dist")  # Type of dependance between variables
  XSpars <- list(list(), list(), list())  # Parameters for dependance
  
  Ysigma <- .2    # Noise amplitude
})
attach(parameters)

#---- Derivative parameters and verifications
p <- length(Alpha)
pvec <- sapply(Alpha, length)
d <- sum(pvec)
pind <- rep(1:p, pvec)

#---- Data generation
# X variables
#XSpars <- Map(c, XSpars, p = pvec)
#Sigma  <- Map(do.call, XSfun, XSpars) # Construction of covariance matrices
#mu <- lapply(pvec, rep, x = 0)
#X <- Map(mvrnorm, n = n, mu = mu, Sigma = Sigma)
#Xall <- Reduce(cbind, X)
X <- sapply(pvec * n, runif)
X <- Map(matrix, X, ncol = pvec)
Xall <- Reduce(cbind, X)

# Indices
AlphaNor <- lapply(Alpha, normalize)
Z <- mapply("%*%", X, AlphaNor)

# Ridge functions
Gpars <- Map(c, Gpars, z = apply(Z, 2, list))
G <- mapply(do.call, Gfuns, Gpars)
G <- scale(G)

# True Y
Y <- Beta0 + G %*% Beta1 

# Simulate Y with noise
Ysim <- replicate(ns, Y + Ysigma * rnorm(n), simplify = F)

#-------------------------------------------
#            Apply models
#-------------------------------------------

# Initialize cluster for parallel computation
cl <- makeCluster(10)

# Load packages in clusters
clusterEvalQ(cl, {
  library(MASS)
  library(mgcv)
  library(scam)
  library(RColorBrewer)
  library(fda)
  library(forecast)
  library(tictoc)
  library(gratia)
  library(quadprog)
  library(quadprogXT)
  library(GA)
  library(scar)
  library(osqp)
  library(np)
  
  source("gaim/Secondary functions.R")
  source("gaim/aim.R")
  source("Benchmark_models.R")
})

# Export useful object
clusterExport(cl, c("Xall", "p", "pind", "pvec", "n"))

# Initialize result objects
mod <- list()
exectime <- list()

#---- AIMs
# Two step
deb <- Sys.time()
mod$TwoStep <- parLapply(cl, Ysim, aim, X, alpha.control = list(norm.type = "sum"), 
  algo.control = list(type = "two.steps"), trace = T)
exectime$TwoStep <- Sys.time() - deb

# Two Step with constraints
deb <- Sys.time()
mod$TwoStepConst <- parLapply(cl, Ysim, aim, X, trace = T,
  algo.control = list(type = "two.steps"), 
  smooth.control = list(shape = c("tp", "mpi", "cx")),
  alpha.control = list(monotone = c(-1, 1, 0), sign.const = c(1, 1, 1),
    norm.type = "sum")
)
exectime$TwoStepConst <- Sys.time() - deb

# Backfitting
deb <- Sys.time()
mod$Backfitting <- parLapply(cl, Ysim, aim, X, 
  alpha.control = list(norm.type = "sum"), 
  algo.control = list(type = "backfitting"), trace = T)
exectime$Backfitting <- Sys.time() - deb

# Backfitting with constraints
deb <- Sys.time()
mod$BackfittingConst <- parLapply(cl, Ysim, aim, X, trace = T, 
  algo.control = list(type = "backfitting"), 
  smooth.control = list(shape = c("tp", "mpi", "cx")),
  alpha.control = list(monotone = c(-1, 1, 0), sign.const = c(1, 1, 1),
    norm.type = "sum")
  )
exectime$BackfittingConst <- Sys.time() - deb

# Gauss-Newton
deb <- Sys.time()
mod$GN <- parLapply(cl, Ysim, aim, X, alpha.control = list(norm.type = "sum"), 
  algo.control = list(type = "gauss.newton"), trace = T)
exectime$GN <- Sys.time() - deb
# GaussNewton with constraints
deb <- Sys.time()
mod$GNConst <- parLapply(cl, Ysim, aim, X, trace = T,
  algo.control = list(type = "gauss.newton"), 
  smooth.control = list(shape = c("tp", "mpi", "cx"), method = "scam"),
  alpha.control = list(monotone = c(-1, 1, 0), sign.const = c(1, 1, 1),
    norm.type = "sum", delta = T)
  )
exectime$GNConst <- Sys.time() - deb

# Fonction optim
deb <- Sys.time()
mod$Optim <- parLapply(cl, Ysim, aim, X, 
  alpha.control = list(norm.type = "sum"), 
  algo.control = list(type = "optim"), trace = T)
exectime$Optim <- Sys.time() - deb
# Fonction optim with constraints
deb <- Sys.time()
mod$OptimConst <- parLapply(cl, Ysim, aim, X, trace = T,
  algo.control = list(type = "optim"), 
  smooth.control = list(shape = c("tp", "mpi", "cx")),
  alpha.control = list(monotone = c(-1, 1, 0), sign.const = c(1, 1, 1),
    norm.type = "sum")
  )
exectime$OptimConst <- Sys.time() - deb

#---- Benchmark models

# PPR
deb <- Sys.time()
mod$PPR <- parLapply(cl, Ysim, function(y){
  res <- ppr(y = y, x = Xall, nterms = p)  # PPR fitting
  alpha <- Map("[", as.data.frame(res$alpha), split(1:9, pind))  # Alphas
  alpha <- Map(normalize, alpha, "sum")
  jf <- 7 + res$smod[1] * (sum(pvec) + 1) # Index for gz
  gz <- matrix(res$smod[jf + 1L:(p * n)], n, p) # gz
  jt <- jf + res$smod[1] * n  # Index for z
  z <- matrix(res$smod[jt + 1L:(p * n)], n, p) # z
  list(alpha = alpha, gz = gz, z = z, coef = c(res$yb, res$beta), 
    fitted = res$fitted.values)
})
exectime$PPR <- Sys.time() - deb

# SCAR (Chen & Samworth, 2015)
deb <- Sys.time()
mod$SCAR <- parLapply(cl, Ysim, function(y){
  res <- scair(x = Xall, y = y, 
    shape = c("l", "in", "cvx"), allnonneg = T)
  gz <- scale(res$componentfit)
  alpha <- Map("[", as.data.frame(res$index), split(1:9, pind))
  alpha <- Map(normalize, alpha, "sum")
  z <- Xall %*% res$index
  coefs <- c(res$constant, attr(gz, "scaled:scale"))
  yhat <- predict(res)
  list(alpha = alpha, gz = gz, z = z, coef = coefs, fitted = yhat)
})
exectime$SCAR <- Sys.time() - deb

# MAVE (Li et al. 2010)
deb <- Sys.time()
mod$MAVE <- parLapply(cl, Ysim, mave, X, 
  alpha.control = list(norm.type = "sum"))
exectime$MAVE <- Sys.time() - deb

# Xia et al. (2005)
deb <- Sys.time()
mod$Xia <- parLapply(cl, Ysim, backfit_llr, X, 
  alpha.control = list(norm.type = "sum", monotone = c(-1, 1, 0)),
  smooth.control = list(shape = c("", "mpi", "cx")))
exectime$Xia <- Sys.time() - deb

#---- Save Results
save(mod, exectime, Ysim, X, Y, G, parameters, 
  file = "Simulation_results.RData")

