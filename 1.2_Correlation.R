###############################################################################
#
#                          Simulation study
#                       Correlation variation
#
###############################################################################

library(doParallel)
library(MASS)
library(Matrix)
library(sfsmisc)
library(abind)

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

# Sample size
n <- 1000

# Number of simulations
ns <- 1000

# Tested correlations
rhovec <- seq(0, 0.75, by = 0.25)

#----- Derived objects -----
nr <- length(rhovec)
p <- length(Alpha)
pvec <- sapply(Alpha, length)
ptot <- sum(pvec)

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

# Save results
results <- vector("list", nr)

# Loop over simulation designs
for (k in 1:nr){
  # Trace computation
  cat("\n", "scenario = ", k, "\n",
    file = "temp/logsim2.txt", append = T)
  
  set.seed(222 + k)
  
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
  
  #----- Apply algorithms
  results[[k]] <- foreach(y = iter(Ysim), i = iter(seq_along(Ysim)),
    .packages = c("cgaim", "Matrix", "scar", "quadprog", "splines2", "scam"), 
    .combine = abind) %dopar% 
  {
    # Initialize data and result object
    y <- Ysim[[i]]
    dat <- c(list(y = y), X)
    alpha_est <- list()
    
    # Apply GAIM
    res <- cgaim(y ~ g(X1) + g(X2) + g(X3), data = dat)
    alpha_est$GAIM <- abs(unlist(res$alpha))
    
    # Apply CGAIM
    res <- cgaim(y ~ 
        g(X1, fcons = "inc", acons = list(monotone = -1, sign.const = 1)) + 
        g(X2, fcons = "inc", acons = list(monotone = 1, sign.const = 1)) + 
        g(X3, fcons = "cvx", acons = list(sign.const = 1)),
      data = dat, smooth_control = list(sp = rep(0, p)))
    alpha_est$CGAIM <- unlist(res$alpha)
    
    # Apply misspecified GAIM
    res <- cgaim(y ~ 
        g(X1, fcons = "inc", acons = list(monotone = 1, sign.const = 1)) + 
        g(X2, fcons = "inc", acons = list(monotone = -1, sign.const = 1)) + 
        g(X3, fcons = "ccv", acons = list(sign.const = 1)),
      data = dat, smooth_control = list(sp = rep(0, p)))
    alpha_est$MGAIM <- unlist(res$alpha)
    
    # Apply cumulative effect GAM (Xia & Tong, 2005)
    x <- X; x[[1]] <- x[[1]][,pvec[1]:1]
    # Rprof(tmp <- tempfile())
    res <- facts(x, y, gshape = c("inc", "inc", "cvx"), 
      thetashape = c("inc", "inc", "nc"))
    # Rprof()
    # summaryRprof(tmp)
    # unlink(tmp)
    alpha <- res$theta; alpha[[1]] <- rev(alpha[[1]])
    alpha_est$FACTS <- unlist(alpha)
    
    # Apply SCAIR
    # res <- scair_apply(y, X, shape = c("in", "in", "cvx"), 
    #   norm.type = "1", epsilon = 1e-3, iter = 50)
    # alpha_est$SCAIR <- abs(res$alpha)
    
    # Apply gMAVE
    res <- gmave(y, X, alpha.control = list(norm.type = "1"))
    alpha_est$gMAVE <- abs(res)
    
    # Apply PPR
    res <- ppr_apply(X, y, norm.type = "1")
    alpha_est$PPR <- abs(res)
    
    # Trace iteration 
    if(i %% 20 == 0) cat("\n", "iter = ", i, as.character(Sys.time()), "\n",
      file = "temp/logsim2.txt", append = T)
    
    # Reorganize into a matrix and return
    array(unlist(alpha_est), dim = c(ptot, length(alpha_est), 1), 
      dimnames = list(NULL, names(alpha_est), NULL))
  }
}

stopCluster(cl)

#---- Save Results
save(results, file = "Results/1.2_Correlation.RData")

# results[2:4] <- lapply(results[2:4], "[",,c(2:1,3:6),)
# results <- c(results2[5], results)

#-------------------------------------------
#    Result summary
#-------------------------------------------

# Method list
method_list <- dimnames(results[[1]])[[2]]
nm <- length(method_list)

# Palette and symbols
pal <- c(rev(brewer.pal(3, "Blues")), 
  colorRampPalette(c("saddlebrown", "white"))(4)[-4])
pchs <- rep(16:17, each = 3)

#----- RMSE in alpha estimation

# Compute RMSE for each realisation
alpha_rmse <- lapply(results, apply, 2:3, function(x){
  sqrt(mean((unlist(x) - unlist(Alpha))^2))
})

# Compute the overall RMSE
rmse_all <- sapply(alpha_rmse, apply, 1, mean)

# Compute RMSE se
rmse_se <- sapply(alpha_rmse, apply, 1, sd)# / sqrt(ns)

# Put everything in data.frame
rmse_df <- data.frame(
  method = factor(rep(method_list, nr), levels = method_list),
  rho = rep(rhovec, each = nm),
  rmse = c(rmse_all), se = c(rmse_se))

# #----- Compare RMSE for null Alphas and non-null Alphas
# 
# # Null Alphas
# zero_rmse <- lapply(results, apply, 2:3, function(x){
#   sqrt(mean((unlist(x) - unlist(Alpha))[unlist(Alpha) == 0]^2))
# })
# rmse_df$rmse_zero <- c(sapply(zero_rmse, apply, 1, mean))
# rmse_df$se_zero <- c(sapply(zero_rmse, apply, 1, sd) / sqrt(ns))
# 
# # NonNull Alphas
# nonzero_rmse <- lapply(results, apply, 2:3, function(x){
#   sqrt(mean((unlist(x) - unlist(Alpha))[unlist(Alpha) != 0]^2))
# })
# rmse_df$rmse_nonzero <- c(sapply(nonzero_rmse, apply, 1, mean))
# rmse_df$se_nonzero <- c(sapply(nonzero_rmse, apply, 1, sd) / sqrt(ns))

#----- Plots

# Create common plot layout
plotlay <- ggplot(rmse_df, aes(x = rho, group = method, 
  col = method, fill = method, shape = method)) + 
  theme_classic() + 
  geom_ribbon(alpha = .2, col = NA) +
  geom_line(size = 1) + 
  geom_point(size = 3) +
  scale_fill_manual(values = pal, name = "Method") + 
  scale_color_manual(values = pal, name = "Method") + 
  scale_shape_manual(values = pchs, name = "Method") +
  scale_y_continuous(name = "RMSE", #trans = "log10", 
    limits = c(0, .3)) + 
  scale_x_continuous(name = "Correlation between variables",
    breaks = rhovec) + 
  theme(panel.grid.major.y = element_line(colour = "grey", linetype = 2))

# Plot overall
plotlay + aes(y = rmse, ymin = rmse - se, ymax = rmse + se)
# plotlay + aes(y = rmse, ymin = rmse - 1.96 * se, ymax = rmse + 1.96 * se)

# # Only nonzero alphas
# plotlay + aes(y = rmse_zero, ymin = rmse_zero - 1.96 * se_zero, 
#   ymax = rmse_zero + 1.96 * se_zero)
# 
# # Only nonzero alphas
# plotlay + aes(y = rmse_nonzero, ymin = rmse_nonzero - 1.96 * se_nonzero, 
#   ymax = rmse_nonzero + 1.96 * se_nonzero)

ggsave("Figures/Figure1b.pdf")

# #----- One example of Alphas
# 
# # Which to show
# shown <- 0.5
# 
# # Restructure into a df
# ind <- which(rhovec == shown)
# alphadf <- expand.grid(alpha = 1:ptot, method = method_list, real = 1:ns)
# alphadf$index <- rep(1:p, pvec)[alphadf$alpha]
# alphadf$alpha <- alphadf$alpha - c(0, cumsum(pvec))[alphadf$index]
# alphadf$val <- c(results[[ind]])
# 
# # True alphas
# truedf <- data.frame(index = rep(1:p, pvec), alpha = c(sapply(pvec, seq_len)),
#   val = unlist(Alpha))
# 
# # Plot
# ggplot(alphadf) + 
#   theme_classic() + 
#   theme(panel.grid.major.y = element_line(colour = grey(.95)),
#     axis.title.x = element_blank(),
#     axis.text.x = element_blank(),
#     axis.ticks.x = element_blank(),
#     axis.line.x = element_blank(),
#     panel.border = element_rect(colour = grey(.9), fill = NA)
#     # panel.background = element_rect(fill = grey(.98))
#   ) + 
#   geom_hline(data = truedf, aes(yintercept = val), size = 1) +
#   geom_boxplot(aes(x = method, y = val, col = method)) + 
#   # geom_hline(yintercept = 0) +
#   facet_grid(index ~ alpha, 
#     labeller = label_bquote(
#       rows = bold(alpha)[.(index)], 
#       cols = alpha[j * .(alpha)])) + 
#   scale_color_manual(values = pal) + 
#   scale_y_continuous(limits = c(0, 1), name = "Estimations")