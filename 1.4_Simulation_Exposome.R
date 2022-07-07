################################################################################
#
#  R code for the simulation study of 
#
#   Masselot et al., 2022
#   Constrained groupwise additive index models
#   Biostatistics
#
#   Section 4.4 Exposome
#
#   Author: Pierre Masselot
#
################################################################################

#-------------------------------------------
# Packages
#-------------------------------------------

#----- Used packages
library(doParallel) # Parallel computing
library(MASS) # For multivariate normal generation
library(ggplot2) # Plotting
library(patchwork) # Assemble ggplots
library(RColorBrewer) # Colorpalette Blues

#----- cgaim package
library(cgaim)

#----- Load correlation matrix (Robinson et al. 2018)
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

set.seed(4)

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

#----- Generate alphas
alphasim <- lapply(nnn, function(nn) replicate(ns, {
  # Draw one non-null variable per group
  indsim <- sapply(split(1:ptot, grpinds), sample, size = 1)
  indexcl <- sapply(split((1:ptot)[-indsim], grpinds[-indsim]), 
    sample, size = 1)
  
  # Draw the rest of non-null variables
  indsim <- c(indsim, sample((1:ptot)[-c(indsim, indexcl)], nn - p))
  
  # Draw alphas: sign with half probability
  alphas <- rep(0, ptot)
  alphas[indsim] <- sample(c(-1, 1), nn, replace = T)
  alphas
}))

#----- Generate outcomes
Ysim <- lapply(seq_along(alphasim), function(a){
  # Generate linear predictor
  linpred <- Xall %*% alphasim[[a]]
  
  # Determine noise for each iteration
  r2a <- nnn[a] * r2
  sig2 <- apply(linpred, 2, var) * r2a / (100 - r2a)
  
  # Generate noise
  epsilon <- mvrnorm(n, rep(0, ns), diag(sig2))
  
  # Final generated response
  linpred + epsilon
})

#-------------------------------------------
#     Simulations
#-------------------------------------------

# Initialize cluster for parallel computation
cl <- makeCluster(max(1, detectCores() - 2))
registerDoParallel(cl)

# To trace simulations
writeLines(c(""), "temp/logsim4.txt")
cat(as.character(as.POSIXct(Sys.time())), file = "temp/logsim4.txt", 
  append = T)

# Packages
packs <- c("cgaim")

#----- Apply CGAIM on each simulation
results <- foreach(k = seq_along(nnn), .packages = packs, .combine = rbind) %:% 
  foreach(i = seq_len(ns), .packages = packs, .combine = rbind) %dopar% 
{
  
  if(i %% 20 == 0) cat("\n", "scenario =", k, "iter = ", i, 
    as.character(Sys.time()), "\n", file = "temp/logsim4.txt", append = T)
  
  # Data preparation
  dat <- c(list(y = Ysim[[k]][,i]), X)
  alpha_est <- list()
  
  # Number of positive alphas in each group
  npos <- tapply(alphasim[[k]][,i], grpinds, function(x) sum(x != 0))
  
  # Apply GAIM
  formula <- sprintf("y ~ %s", 
    paste(sprintf("g(%s)", 
      names(groups)), collapse = " + "))
  res <- cgaim(as.formula(formula), data = dat,
    control = list(sm_pars = list(sp = rep(0, p))))
  alpha_est$GAIM <- unlist(mapply("*", res$alpha, npos))
  
  # Apply CGAIM
  formula <- sprintf("y ~ %s", paste(sprintf("g(%s, fcons = '%s')", 
    names(groups), c("cvx", "inc", "inc", "inc", "inc")), 
    collapse = " + "))
  res <- cgaim(as.formula(formula), data = dat, 
    Cmat = rbind(diag(ptot), -diag(ptot)), bvec = rep(rep(-1/npos, pvec), 2), 
    control = list(sm_pars = list(sp = rep(0, p))))
  
  # Correct for numerical errors in OSQP/quadprog
  alpha_est$CGAIM <- pmin(pmax(unlist(mapply("*", res$alpha, npos)), -1), 1)
  
  # Reorganize into a data.frame and return
  data.frame(nnn = nnn[k], simu = i, model = rep(names(alpha_est), each = ptot),
    alpha = rep(1:ptot, length(alpha_est)), est = unlist(alpha_est), 
    true = rep(alphasim[[k]][,i], length(alpha_est)))
}

# Stop parallel
stopCluster(cl)

#-------------------------------------------
#    Result summary
#-------------------------------------------

# Models applied
method_list <- unique(results$model)
nm <- length(method_list)

# Compute mean and confidence intervals
ovmean <- aggregate(est ~ model + true + nnn, results, mean)
names(ovmean)[4] <- "mean"
ovlow <- aggregate(est ~ model + true + nnn, results, quantile, .025)
names(ovlow)[4] <- "low"
ovhigh <- aggregate(est ~ model + true + nnn, results, quantile, .975)
names(ovhigh)[4] <- "high"

# Put together in data.frame
ovres <- Reduce(merge, list(ovmean, ovlow, ovhigh))

#----- Figure 4: Plot alpha distribution

# Method list
method_list <- unique(results$model)
nm <- length(method_list)

# Palette for models
pal <- tail(brewer.pal(pmax(nm, 3), "Blues"), nm)
names(pal) <- method_list

# Plot
ggplot(ovres, aes(x = nnn, group = interaction(true, model), col = model)) + 
  theme_classic() + 
  geom_hline(yintercept = c(-1, 0, 1), linetype = 2) + 
  geom_errorbar(aes(ymin = low, ymax = high), width = 1, size = .5,
    position = position_dodge(width = 1)) +
  geom_point(aes(y = mean), position = position_dodge(width = 1), 
    size = 3) +
  scale_x_continuous(name = "Number of true predictors", breaks = nnn) + 
  scale_y_continuous(name = "Estimated", n.breaks = 6) + 
  scale_color_manual(name = "", values = pal)
