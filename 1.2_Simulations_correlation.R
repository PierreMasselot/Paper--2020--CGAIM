###############################################################################
#
#                          Simulation study
#                       Sample size variation
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
source("1.0_Simulation_functions.R")
source("1.0_Parameters.R")

#-------------------------------------------
#     Simulations
#-------------------------------------------

# Tested correlations
corvec <- seq(0, .8, by = .2)
nc <- length(corvec)

# Initialize cluster for parallel computation
cl <- makeCluster(2)
# Transfer objects in cluster
clusterExport(cl, ls())
clusterEvalQ(cl, {
  library(parallel)
  library(MASS)
  library(Matrix)
  library(cgaim)
})

# Initialize result objects
dat_cor <- results <- vector("list", nc) # Generated data

# Save results
time_cor <- rep(list(vector("numeric", nc)), nmod)
g_cor <- rep(list(vector("list", nc)), nmod)
z_cor <- rep(list(vector("list", nc)), nmod)
alpha_cor <- rep(list(vector("list", nc)), nmod)
yhat_cor <- rep(list(vector("list", nc)), nmod)

# Loop over simulation designs
for (k in 1:nc){
  print(k); flush.console()
  
  #----- Generate data  
  dat_cor[[k]] <- generate_data(n, ns, Alpha, Gfuns, Gpars,
    Beta0, Beta1, corvec[k], Ysigma)
  Xall <- Reduce(cbind, dat_cor[[k]]$X) # Useful for PPR
  
  # Transfer objects in cluster
  clusterExport(cl, c("Xall", "dat_cor", "k", "corvec"))
  
  #---- CGAIM models
  ## Unconstrained GAIM
  print("GAIM"); flush.console()
  deb <- Sys.time()
  results <- parLapply(cl, dat_cor[[k]]$Y, function(y){
    dat <- c(list(y = y), dat_cor[[k]]$X)
    cgaim(y ~ g(X1) + g(X2) + g(X3), data = dat, 
      alpha.control = list(norm.type = "sum"))
  })
  time_cor[[1]][k] <- Sys.time() - deb
  # Estimated alphas and functions
  g_cor[[1]][[k]] <- sapply(results, "[[", "gfit", simplify = "array")
  z_cor[[1]][[k]] <- sapply(results, "[[", "indexfit", simplify = "array")
  alphas <- sapply(results, "[[", "alpha")
  alpha_cor[[1]][[k]] <- apply(alphas, 2, unlist)
  yhat_cor[[1]][[k]] <- sapply(results, "[[", "fitted")
  
  ## CGAIM
  print("CGAIM"); flush.console()
  deb <- Sys.time()
  # Apply model
  results <- parLapply(cl, dat_cor[[k]]$Y, function(y){
    dat <- c(list(y = y), dat_cor[[k]]$X)
    cgaim(y ~ g(X1, bs = "mpi", constraints = list(monotone = -1, sign.const = 1)) + 
        g(X2, bs = "mpi", constraints = list(monotone = 1, sign.const = 1)) + 
        g(X3, bs = "cx", constraints = list(sign.const = 1)),
      data = dat, alpha.control = list(norm.type = "sum"))
  })
  time_cor[[2]][k] <- Sys.time() - deb
  # Estimated alphas and functions
  g_cor[[2]][[k]] <- sapply(results, "[[", "gfit", simplify = "array")
  z_cor[[2]][[k]] <- sapply(results, "[[", "indexfit", simplify = "array")
  alphas <- sapply(results, "[[", "alpha")
  alpha_cor[[2]][[k]] <- apply(alphas, 2, unlist)
  yhat_cor[[2]][[k]] <- sapply(results, "[[", "fitted")
  
  ## Misspecified cGAIM
  print("MGAIM"); flush.console()
  deb <- Sys.time()
  results <- parLapply(cl, dat_cor[[k]]$Y, function(y){
    dat <- c(list(y = y), dat_cor[[k]]$X)
    cgaim(y ~ g(X1, bs = "mpi", constraints = list(sign.const = 1)) + 
        g(X2, bs = "mpi", constraints = list(monotone = -1, sign.const = 1)) + 
        g(X3, bs = "cx", constraints = list(monotone = 1, sign.const = 1)),
      data = dat, alpha.control = list(norm.type = "sum"))
  })
  time_cor[[3]][k] <- Sys.time() - deb
  # Estimated alphas and functions
  g_cor[[3]][[k]] <- sapply(results, "[[", "gfit", simplify = "array")
  z_cor[[3]][[k]] <- sapply(results, "[[", "indexfit", simplify = "array")
  alphas <- sapply(results, "[[", "alpha")
  alpha_cor[[3]][[k]] <- apply(alphas, 2, unlist)
  yhat_cor[[3]][[k]] <- sapply(results, "[[", "fitted")  
  
  #---- Benchmark models
  ## PPR
  print("PPR"); flush.console()
  deb <- Sys.time()
  results <- parLapply(cl, dat_cor[[k]]$Y, function(y){
    res <- ppr(y = y, x = Xall, nterms = p)  # PPR fitting
    alpha <- Map("[", as.data.frame(res$alpha), split(1:9, pind))  # Alphas
    alpha <- Map(cgaim:::normalize, alpha, "sum")
    n <- corvec[k]
    jf <- 7 + res$smod[1] * (sum(pvec) + 1) # Index for gz
    gz <- matrix(res$smod[jf + 1L:(p * n)], n, p) # gz
    jt <- jf + res$smod[1] * n  # Index for z
    z <- matrix(res$smod[jt + 1L:(p * n)], n, p) # z
    list(alpha = alpha, gz = gz, z = z, coef = c(res$yb, res$beta), 
      fitted = res$fitted.values)
  })
  time_cor[[4]][k] <- Sys.time() - deb
  # Estimated alphas and functions
  g_cor[[4]][[k]] <- sapply(results, "[[", "gz", simplify = "array")
  z_cor[[4]][[k]] <- sapply(results, "[[", "z", simplify = "array")
  alphas <- sapply(results, "[[", "alpha")
  alpha_cor[[4]][[k]] <- apply(alphas, 2, unlist)
  yhat_cor[[4]][[k]] <- sapply(results, "[[", "fitted")
  
  
  ## gMAVE
  print("gMAVE"); flush.console()
  deb <- Sys.time()
  results <- parLapply(cl, dat_cor[[k]]$Y, mave, dat_cor[[k]]$X, 
    alpha.control = list(norm.type = "sum"))
  time_cor[[5]][k] <- Sys.time() - deb
  # Estimated alphas and functions
  g_cor[[5]][[k]] <- sapply(results, "[[", "gz", simplify = "array")
  z_cor[[5]][[k]] <- sapply(results, "[[", "z", simplify = "array")
  alphas <- sapply(results, "[[", "alpha")
  alpha_cor[[5]][[k]] <- apply(alphas, 2, unlist)
  yhat_cor[[5]][[k]] <- sapply(results, "[[", "fitted")
}

stopCluster(cl)

#---- Save Results
save(corvec, dat_cor, time_cor, g_cor, alpha_cor, z_cor, yhat_cor, 
  file = "Results/1.2_Simulations_correlation.RData")


#-------------------------------------------
#     Plots
#-------------------------------------------

load("Results/1.2_Simulations_correlation.RData")

### MSE and MISE plot

# Functions error (MISE)
trueGs <- lapply(dat_cor, "[[", "G")
trueZs <- lapply(dat_cor, function(x){
  mapply("%*%", x$X, Alpha)
})

ISEs <- array(NA, dim = c(nc, p, nmod, ns))
for (e in 1:nc){
  for (j in 1:p){
    trueFunc <- splinefun(trueZs[[e]][,j], trueGs[[e]][,j])
    trueEval <- trueFunc(seq(min(trueZs[[e]][,j]), max(trueZs[[e]][,j]), 
      length.out = 1000))
    for (m in 1:nmod){
      modFunc <- Map(function(z, g){
          if (all(is.na(g))){
            out <- rep(NA, 1000)
          } else{
            out <- splinefun(z, g)(seq(min(z), max(z), length.out = 1000))
          }
          return(out)
        }, 
        as.data.frame(z_cor[[m]][[e]][,j,]),
        as.data.frame(g_cor[[m]][[e]][,j,]))
      # Integrate the difference by trapezoid approximation
      ISEs[e,j,m,] <- mapply(function(z, g){
          ifelse(all(is.na(g)) || all(g == 0), NA, 
            integrate.xy(seq(min(z), max(z), length.out = 1000), 
              (g - trueEval)^2))
        }, as.data.frame(z_cor[[m]][[e]][,j,]), modFunc) 
    }
  }
}
MISEs <- apply(ISEs, 1:3, mean, na.rm = T)


# Alpha errors
true_alphas <- unlist(Alpha)

# Compute RMSEs for each alpha
alpha_errors <- lapply(alpha_cor, lapply, apply, 2, "-", true_alphas)
alpha_mses <- lapply(alpha_errors, sapply, apply, 1, 
  function(x) sqrt(sum(x^2, na.rm = T)))

alpha_mse_ses <- Map(function(er, mse){
  mapply(function(er1, mse1){
    vars <- apply(er1, 2, function(x) x^2 - mse1)
    rowSums(vars) / (ns * (ns - 1))
  }, er, as.data.frame(mse))
}, alpha_errors, alpha_mses)

alpha_mses_byZ <- lapply(alpha_mses, aggregate, by = list(Group = pind), mean)
error_rg <- range(sapply(alpha_mses, function(x) range(x[,-1])))

# plot
x11(height = 10, width = 15)
layout(rbind(matrix(1:(2*p), nrow = 2), (2*p) + 1), 
  heights = c(.45, .45, .1))
for (j in 1:p){
  group_errors <- sapply(alpha_mses_byZ, "[", j, -1)
  matplot(corvec, group_errors, type = "b", lwd = 3,
    col = mod_pal, pch = mod_pch, log = "y", xlab = "Correlation coefficient",
    main = bquote(.(letters[j]) * ") Group" ~ .(j) * ":" ~ italic(alpha)),
      #paste0(letters[j], ") Group ", j, expression(alpha)),
    ylim = error_rg, 
    ylab = ifelse(j == 1, "MSE", ""), cex.lab = 1.5, cex.axis = 1.3, 
    cex.main = 2, cex = 1.5, lty = mod_lty)
  matplot(corvec, MISEs[,j,], type = "b", lwd = 3,
    col = mod_pal, pch = mod_pch, log = "y", xlab = "Correlation coefficient",
    main = bquote(.(letters[p + j]) * ") Group" ~ .(j) * ":" ~ italic(g)), 
    ylim = range(MISEs), 
    ylab = ifelse(j == 1, "MISE", ""), cex.lab = 1.5, cex.axis = 1.3, 
    cex.main = 2, cex = 1.5, lty = mod_lty)
}
par(mar = rep(0,4))
plot.new()
legend("center", mod_names, col = mod_pal, lwd = 2, lty = mod_lty, 
  ncol = nmod, cex = 1.5, pch = mod_pch)
  
dev.print(png, filename = "Results/1.2_Fig3.png", res = 100, 
  width = dev.size()[1], height = dev.size()[2], units = "in")
dev.copy2eps(file = "Results/1.2_Fig3.eps")

