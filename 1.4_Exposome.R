###############################################################################
#
#                          Simulation study
#                       Sample size variation
#
###############################################################################

library(foreach)
library(doParallel)
library(MASS)
library(Matrix)
library(sfsmisc)
library(scar)
library(ggplot2)
library(abind)
library(RColorBrewer)


# Should be installed from github
# install_github("PierreMasselot/cgaim")
library(cgaim)
source("0_Useful_functions.R")
source("1.0_Benchmark_models.R")

# Load correlation matrix (Robinson et al. 2018)
corrmat <- data.matrix(read.table("Data/1.4_ExposomeCorrelation.csv", 
  sep = ","))

#-------------------------------------------
#     Parameters
#-------------------------------------------

# Number of observations
n <- 1200

# Number of simulations
ns <- 1000

# Group specification
groups <- list(meteorology = c(1:4), ap = c(5:9), road = c(10, 16:18), 
  natural = c(11:15), built = c(19:28))

# Number of non-null elements
nnn <- c(5, 10, 15)

# R2 baseline (%)
r2 <- 3

#-------------------------------------------
#     Exposome generation
#-------------------------------------------

# Reorder correlation matrix according to groups
corrmat <- corrmat[unlist(groups), unlist(groups)]

# Number of variables and groups
p <- length(groups)
pvec <- lengths(groups)
ptot <- ncol(corrmat)
grpinds <- rep(1:p, pvec)

# Generate Exposome from correlation matrix
Xall <- mvrnorm(n, rep(0, ptot), corrmat)
colnames(Xall) <- colnames(corrmat)

# Create list by groups
X <- tapply(1:ptot, grpinds, function(i) Xall[,i])
names(X) <- names(groups)

#-------------------------------------------
#     Simulations
#-------------------------------------------

# Initialize cluster for parallel computation
cl <- makeCluster(max(1, detectCores() - 2))
registerDoParallel(cl)

# To trace simulations
writeLines(c(""), "temp/logsim5.txt")
cat(as.character(as.POSIXct(Sys.time())), file = "temp/logsim5.txt", 
  append = T)

results <- alphasim <- vector("list", length(nnn))

# Loop over simulation designs
for (k in seq_along(nnn)){
  
  cat("\n", "scenario = ", k, "\n",
    file = "temp/logsim5.txt", append = T)
  
  set.seed(555 + k)
  
  #----- Generate response  
  
  # Generate non-null alphas for each iteration
  alphasim[[k]] <- replicate(ns, {
    indsim <- sapply(split(1:ptot, grpinds), sample, size = 1)
    indsim <- c(indsim, sample((1:ptot)[-indsim], nnn[k] - p))
    alphas <- rep(0, ptot)
    alphas[indsim] <- 1
    alphas
  })
  
  # Generate linear predictor
  linpred <- Xall %*% alphasim[[k]]
  
  # Determine noise for each iteration
  sig2 <- apply(linpred, 2, var) * r2 / (100 - r2)
  
  # Generate noise
  epsilon <- mvrnorm(n, rep(0, ns), diag(sig2))
  
  # Final generated response
  Ysim <- linpred + epsilon
  
  #----- Apply CGAIM on each simulation
  results[[k]] <- foreach(i = iter(seq_len(ns)),
    .packages = c("cgaim", "Matrix", "glmnet"), 
    .combine = abind) %dopar% 
  {
    
    # Data preparation
    dat <- c(list(y = Ysim[,i]), X)
    alpha_est <- list()
    
    # Number of positive alphas in each group
    npos <- tapply(alphasim[[k]][,i], grpinds, sum)
    
    # Apply GAIM
    formula <- sprintf("y ~ %s", 
      paste(sprintf("g(%s)", 
        names(groups)), collapse = " + "))
    res <- cgaim(as.formula(formula), data = dat,
      smooth_control = list(sp = rep(0, p)))
    alpha_est$GAIM <- unlist(mapply("*", res$alpha, npos))
    
    # Apply CGAIM
    formula <- sprintf("y ~ %s", paste(sprintf("g(%s, fcons = '%s')", 
      names(groups), c("cvx", "inc", "inc", "inc", "inc")), 
      collapse = " + "))
    res <- cgaim(as.formula(formula), data = dat, 
      alpha_control = list(Cmat = diag(ptot)), 
      smooth_control = list(sp = rep(0, p)))
    alpha_est$CGAIM <- unlist(mapply("*", res$alpha, npos))
    
    # Apply Elastic net
    # res <- cv.glmnet(Xall, Ysim[,i], relax = T)
    # alpha_est$ENET <- as.matrix(coef(res, s = "lambda.1se")[-1])
    
    # Apply GMAVE
    # res <- gmave(Ysim[,i], X, alpha.control = list(norm.type = "L1"))
    # alpha_est$GMAVE <- mapply("*", split(res, rep(1:p, pvec)), npos)
    
    # Trace iteration 
    if(i %% 20 == 0) cat("\n", "iter = ", i, as.character(Sys.time()), "\n",
      file = "temp/logsim5.txt", append = T)
    
    # Reorganize into a matrix and return
    array(unlist(alpha_est), dim = c(ptot, length(alpha_est), 1), 
      dimnames = list(NULL, names(alpha_est), NULL))
  }
}

stopCluster(cl)

#---- Save Results
save(results, alphasim,
  file = "Results/1.5_Exposome.RData")

#-------------------------------------------
#    Result summary
#-------------------------------------------

# Models applied
method_list <- dimnames(results[[1]])[[2]]
nm <- length(method_list)

# Data.frame with estimated alphas
alphadf <- expand.grid(var = 1:ptot, sim = 1:ns, method = method_list, 
  nonnull = nnn)
alphadf$alpha <- c(sapply(results, function(x) c(aperm(x, c(1, 3:2)))))
alphadf$true <- c(sapply(alphasim, rep, nm))

#----- Plot all alphas

# Compute mean and confidence intervals
ovmean <- aggregate(alpha ~ method + true + nonnull, alphadf, mean)
names(ovmean)[4] <- "mean"
ovlow <- aggregate(alpha ~ method + true + nonnull, alphadf, quantile, .025)
names(ovlow)[4] <- "low"
ovhigh <- aggregate(alpha ~ method + true + nonnull, alphadf, quantile, .975)
names(ovhigh)[4] <- "high"

# Put together in data.frame
ovres <- Reduce(merge, list(ovmean, ovlow, ovhigh))

# Plot
ggplot(ovres, aes(x = nonnull, group = method, col = method)) + 
  theme_classic() + 
  geom_pointrange(aes(y = mean, ymin = low, ymax = high), 
    position = position_dodge(width = 2), size = .9) +
  scale_x_continuous(name = "Number of true predictors", breaks = nnn) + 
  scale_y_continuous(name = "Estimated", n.breaks = 6) + 
  scale_color_manual(name = "", values = brewer.pal(3, "Blues")[-1]) + 
  geom_hline(yintercept = c(0, 1), linetype = 2)

ggsave("Figures/Figure3.pdf")
