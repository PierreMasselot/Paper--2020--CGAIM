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
  function(z, lambda = 1) log(z - min(z) + .1), 
  function(z, lambda = 5) 1 / (1 + exp(-lambda * scale(z))), 
  function(z, a1 = 50, a2 = .8, b1 = .5, b2 = 1.5) 
    a1 * exp(-b1 * scale(z)) + a2 * exp(b2 * scale(z))
)

# Number of simulations
ns <- 1000

# Number of observations
n <- 1000

# Number of nonull indices
pp <- 1:2

# Indices specification for CGAIM
models <- list("g(X1, fcons = 'inc', acons = list(monotone = -1, sign.const = 1))",
  "g(X2, fcons = 'inc', acons = list(monotone = 1, sign.const = 1))",
  "g(X3, fcons = 'cvx', acons = list(sign.const = 1))"
)

#----- Derived objects -----
p <- length(Alpha)
pvec <- sapply(Alpha, length)
ptot <- sum(pvec)
np <- length(pp)

#-------------------------------------------
#     Generate indices
#-------------------------------------------

set.seed(333)

# Uncorrelated variables
X <- Map(MASS::mvrnorm, n = n, mu = lapply(pvec, rep, x = 0), 
  Sigma = lapply(pvec, diag))
names(X) <- sprintf("X%i", 1:p)

# Create indices from generated variables
Z <- Map("%*%", X, Alpha)

# Apply ridge functions
G <- mapply(function(g, z) do.call(g, list(z = z)), Gfuns, Z)

#-------------------------------------------
#     Simulations
#-------------------------------------------

# Initialize cluster for parallel computation
cl <- makeCluster(max(1, detectCores() - 2))
registerDoParallel(cl)

# To trace simulations
writeLines(c(""), "temp/logsim3.txt")
cat(as.character(as.POSIXct(Sys.time())), file = "temp/logsim3.txt", 
  append = T)

# Store results
results <- nonnull <- vector("list", np)

# Loop over simulation designs
for (k in 1:np){
  
  cat("\n", "scenario = ", k, "\n",
    file = "temp/logsim3.txt", append = T)
  
  set.seed(3333 + k)
  
  #----- Generate response
  
  # For each realization choose which are the non null indices
  nonnull[[k]] <- replicate(ns, {
    sim <- rep(0, p)
    sim[sample.int(p, pp[k])] <- 1
    sim
  })
  
  # Generate linear predictor
  Y <- 5 + scale(G) %*% nonnull[[k]]
  
  # Add noise
  Ysim <- Y + matrix(rnorm(n * ns, 0, .2), nrow = n, ncol = ns)
  
  #----- Loop on simulations
  results[[k]] <- foreach(i = iter(seq_len(ns)),
    .packages = c("cgaim", "Matrix", "scar", "quadprog", "splines2", "scam"), 
    .combine = abind) %dopar% 
  {
    # Initialize data and result object
    dat <- c(list(y = Ysim[,i]), X)
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

    # Trace iteration 
    if(i %% 20 == 0) cat("\n", "iter = ", i, as.character(Sys.time()), "\n",
      file = "temp/logsim3.txt", append = T)
    
    # Reorganize in array and return
    array(selectvar, dim = c(dim(selectvar), 1), 
      dimnames = c(dimnames(selectvar), list(NULL)))
  }
}

stopCluster(cl)

#---- Save Results
save(results, nonnull, file = "Results/1.3_IndexSelection.RData")

#-------------------------------------------
#    Result summary
#-------------------------------------------

# Method list
method_list <- dimnames(results[[1]])[[1]]
nm <- length(method_list)

# Palette and symbols
pal <- c(brewer.pal(3, "Blues")[-1])

#----- Compute detection scores

# Compute sensitivity for each realisation
sensitivity <- mapply(function(res, nn, p) 
  apply(res, 1, function(x)  colSums(x & nn)) / p,
  results, nonnull, pp)

# Compute specificity for each realisation
specificity <- mapply(function(res, nn, p) 
  apply(res, 1, function(x)  colSums(!x & !nn)) / p,
  results, nonnull, p - pp)

# Put them together in big data.frame
scoredf <- expand.grid(real = 1:ns, method = method_list, ptrue = factor(pp))
scoredf$sensitivity <- c(sensitivity)
scoredf$specificity <- c(specificity)

#----- Ploting them
# 
# # Scatterplot
# 
# ggplot(scoredf) + theme_classic() +
#   geom_point(aes(x = sensitivity, y = specificity, color = method, 
#     shape = ptrue), 
#     size = 3, position = position_jitter(width = .1, height = .1)) +
#   scale_color_manual(name = "Method", values = pal) + 
#   scale_shape_manual(name = "Number of true indices", values = c(16:17)) + 
#   scale_x_continuous(name = "Sensitivity", limits = c(-0.1, 1.1),
#     breaks = c(0, pp) / max(pp)) + 
#   scale_y_continuous(name = "Specificity", limits = c(-0.1, 1.1),
#     breaks = c(0, pp) / max(pp)) + 
#   theme(panel.grid.major = element_line(colour = "grey", linetype = 2))

# Density plot
densitydf <- expand.grid(method = method_list, ptrue = factor(pp),
  sensitivity = c(0, pp) / max(pp), specificity = c(0, pp) / max(pp))
densitydf <- merge(densitydf, aggregate(real ~ ., scoredf, length), all.x = T)  
densitydf$real[is.na(densitydf$real)] <- 0
densitydf$real <- 100 * densitydf$real / ns

labs <- sprintf("%i %s", pp, ifelse(pp == 1, "index", "indices"))
names(labs) <- pp

ggplot(densitydf) + theme_classic() +
  geom_tile(aes(x = sensitivity, y = specificity, fill = real)) +
  scale_fill_continuous(name = "Proportion (%)", limits = c(0, 100),
    low = "white", high = brewer.pal(9, "Blues")[9]) +
  # scale_fill_distiller(name = "Density", limits = c(0, 1), 
  #   palette = "Blues", direction = 1) +
  scale_x_continuous(name = "Sensitivity", #limits = c(-0.5, 1.5),
    breaks = c(0, pp) / max(pp)) + 
  scale_y_continuous(name = "Specificity", #limits = c(-0.5, 1.5),
    breaks = c(0, pp) / max(pp)) + 
  theme(panel.grid.minor = element_line(colour = "grey", linetype = 2), 
    panel.ontop = T, panel.background = element_rect(fill = NA),
    panel.border = element_rect(colour = grey(.9), fill = NA),
    strip.text = element_text(size = 13)) + 
  facet_grid(vars(ptrue), vars(method), labeller = labeller(
    method = label_value, ptrue = labs)) + 
  guides(fill = guide_colourbar(title.vjust = 2))

ggsave("Figures/Figure2.pdf")
