################################################################################
#
#  R code for the simulation study of 
#
#   Masselot et al., 2022
#   Constrained groupwise additive index models
#   Biostatistics
#
#   Section 4.2 Index selection
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
# Should be installed from github
# install_github("PierreMasselot/cgaim")
library(cgaim)

#-------------------------------------------
#     Parameters
#-------------------------------------------

#----- True functions

# True index weights alphas
Alpha <- list(
  c(.7, .2, .1, 0), 
  c(0, 0, .5, .5), 
  c(.2, .4, .3, .1)
)

# True ridge functions G
Gfuns <- c(
  function(z, lambda = 1) log(z - min(z) + .1), 
  function(z, lambda = 5) 1 / (1 + exp(-lambda * scale(z))), 
  function(z, lambda = c(0.2118881, 0.3006585, 0.0982663, 0.0153671, 0.0016265))
    poly(z, 5) %*% lambda
)

#----- Simulation scenarios

# Number of simulations
ns <- 1000

# Number of observations
n <- 1000

# Number of non-null indices
pp <- 1:2

# Noise levels
Sigmas <- c(.2, .5, 1)

# Full grid of scenarios
scenarios <- expand.grid(pp = pp, sigma = Sigmas)

#----- Models

# Index specification for CGAIM
models <- list(
  "g(X1, fcons = 'inc', acons = list(monotone = -1, sign = 1))",
  "g(X2, fcons = 'inc', acons = list(monotone = 1, sign = 1))",
  "g(X3, fcons = 'cvx', acons = list(sign = 1))"
)

#----- Derived objects -----
p <- length(Alpha)
pvec <- sapply(Alpha, length)
ptot <- sum(pvec)
nsce <- nrow(scenarios)

#-------------------------------------------
#     Generate data
#-------------------------------------------

set.seed(2)

#----- Generate linear predictor

# Matrix of uncorrelated predictors
X <- Map(MASS::mvrnorm, n = n, mu = lapply(pvec, rep, x = 0), 
  Sigma = lapply(pvec, diag))
names(X) <- sprintf("X%i", 1:p)

# Create indices from generated variables
Z <- Map("%*%", X, Alpha)

# Apply ridge functions
G <- mapply(function(g, z) do.call(g, list(z = z)), Gfuns, Z)

#----- Generate outcome
datasim <- apply(scenarios, 1, function(s){ 
  
  # Sample which indices are in the model
  nonnull <- replicate(ns, {
    sim <- rep(0, p)
    sim[sample.int(p, s["pp"])] <- 1
    sim
  })
  
  # Generate outcome
  Ysim <- 5 + scale(G) %*% nonnull + matrix(rnorm(n * ns, 0, s["sigma"]), 
    nrow = n, ncol = ns)
  
  # Return
  list(nonnull = nonnull, Ysim = Ysim)
  
}, simplify = F)

#-------------------------------------------
#     Simulations
#-------------------------------------------

# Initialize cluster for parallel computation
cl <- makeCluster(max(1, detectCores() - 2))
registerDoParallel(cl)

# To trace simulations
writeLines(c(""), "temp/logsim2.txt")
cat(as.character(as.POSIXct(Sys.time())), file = "temp/logsim2.txt", 
  append = T)

# Packages
packs <- c("cgaim")

# Loop
results <- foreach(k = seq_len(nsce), .packages = packs, .combine = rbind) %:%
  foreach(i = seq_len(ns), .packages = packs, .combine = rbind) %dopar% 
  {
    
    # Trace iteration 
    if(i %% 20 == 0) cat("\n", "scenario =", k, "iter = ", i, 
      as.character(Sys.time()), "\n", file = "temp/logsim2.txt", append = T)
    
    # Initialize data and result object
    dat <- c(list(y = datasim[[k]]$Ysim[,i]), X)
    selectvar <- matrix(NA, nrow = 0, ncol = p)
    
    #----- GAIM
    selected <- c()
    best <- Inf
    
    # Loop on number of variable for forward steps
    for (v in seq_len(p)){
      remain <- setdiff(seq_len(p), selected)
      gcvs <- rep(NA, length(remain))
      
      # Loop on remaining variables
      for (w in seq_along(remain)){
        indices <- sprintf("g(X%i)", c(selected, remain[w]))
        rhs <- paste(indices, collapse = " + ")
        form <- paste0("y ~ ", rhs)
        res <- cgaim(as.formula(form), data = dat, 
          alpha_control = list(solver = "osqp"))
        gcvs[w] <- res$gcv
      }
      
      # If GCV decreases include variable that leads to largest decrease
      if (any(gcvs < best)){
        sel <- which.min(gcvs)
        best <- gcvs[sel]
        selected <- c(selected, remain[sel])
      } else {
        break
      }
    }
    
    # Store selected variables
    selectvar <- rbind(selectvar, GAIM = 0)
    selectvar["GAIM", selected] <- 1
    
    #----- GAIM
    selected <- c()
    best <- Inf
    
    # Loop on number of variable for forward steps
    for (v in seq_len(p)){
      remain <- setdiff(seq_len(p), selected)
      gcvs <- rep(NA, length(remain))
      
      # Loop on remaining variables
      for (w in seq_along(remain)){
        rhs <- paste(models[c(selected, remain[w])], collapse = " + ")
        form <- paste0("y ~ ", rhs)
        res <- cgaim(as.formula(form), data = dat, 
          alpha_control = list(solver = "osqp"),
          smooth_control = list(sp = rep(0, v)))
        gcvs[w] <- res$gcv
      }
      
      # If GCV decreases include variable that leads to largest decrease
      if (any(gcvs < best)){
        sel <- which.min(gcvs)
        best <- gcvs[sel]
        selected <- c(selected, remain[sel])
      } else {
        break
      }
    }
    
    # Store selected variables
    selectvar <- rbind(selectvar, CGAIM = 0)
    selectvar["CGAIM", selected] <- 1
    
    #----- Return 
    
    # Reorganize in data.frame and return
    data.frame(scenarios[k,,drop = F], simu = i, 
      expand.grid(method = rownames(selectvar), var = 1:p), 
      selected = c(selectvar))
}

# Stop parallel
stopCluster(cl)

#-------------------------------------------
#    Result summary
#-------------------------------------------

# Method list
method_list <- unique(results$method)
nm <- length(method_list)

# Palette for models
pal <- tail(brewer.pal(pmax(nm, 3), "Blues"), nm)
names(pal) <- method_list

# Shapes for true number of variables
shp <- 20 + seq_along(pp)

#----- Compute detection scores

# Add info about true presence
results$true <- rep(c(sapply(datasim, "[[", "nonnull")), each = nm)

# Compute sensitivity
sensitivity <- aggregate(selected ~ pp + sigma + method, 
  data = subset(results, as.logical(true)), mean)
colnames(sensitivity)[4] <- "sensitivity"

# Compute specificity for each realisation
specificity <- aggregate(selected ~ pp + sigma + method, 
  data = subset(results, !as.logical(true)), function(x) 1 - mean(x))
colnames(specificity)[4] <- "specificity"

# Put them together 
scores <- merge(sensitivity, specificity)
scores <- scores[order(scores$sigma, decreasing = T),]

#----- Figure 2: Sensitivity vs specificity

ggplot(scores) + theme_classic() +
  geom_point(aes(x = sensitivity, y = specificity, 
    fill = method, shape = as.factor(pp), size = sigma),
    col = 1) +
  scale_fill_manual(values = pal, labels = c("GAIM", "CGAIM"), 
    name = "Method") + 
  scale_size_area(name = "Error level", breaks = Sigmas, 
    max_size = 8) + 
  scale_shape_manual(values = shp, labels = pp, 
    name = "Number of true indices") + 
  guides(fill = guide_legend("Method", 
      override.aes = list(shape = 21, size = 4)),
    shape = guide_legend("Number of true indices", 
      override.aes = list(size = 4)),
    size = guide_legend("Error level", override.aes = list(fill = "white",
      shape = 21))) +
  lims(x = c(0, 1), y = c(0, 1)) + 
  labs(x = "Sensitivity", y = "Specificity") + 
  theme(legend.position = c(0.1, 0.95), legend.justification = c(0, 1), 
    panel.grid.major = element_line(color = grey(.95)),
    panel.border = element_rect(fill = NA, color = "black"))
  