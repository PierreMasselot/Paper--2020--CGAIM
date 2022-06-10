################################################################################
#
#  R code for the simulation study of 
#
#   Masselot et al., 2022
#   Constrained groupwise additive index models
#   Biostatistics
#
#   Section 4.1 Index estimation
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
library(Matrix) # Used in gMAVE
library(ggplot2) # Plotting
library(patchwork) # Assemble ggplots
library(RColorBrewer) # Colorpalette Blues
library(splines2) # For increasing and convex splines used in facts
library(quadprog) # Quadratic programming used in facts
library(scam) # Utility functions used in facts

#----- cgaim package
# Should be installed from github
# install_github("PierreMasselot/cgaim")
library(cgaim)

#----- Custom functions for benchmark models
source("0.1_gMAVE.R")
source("0.2_FACTS.R")
source("0.3_PPR.R")

#-------------------------------------------
#     Parameters
#-------------------------------------------

#----- True model

# True index weights alphas
Alpha <- list(
  c(.7, .2, .1, 0), 
  c(0, 0, .5, .5), 
  c(.2, .4, .3, .1)
)

# True ridge functions G
Gfuns <- c(
  function(z) log(z - min(z) + .1), 
  function(z, lambda = 20) 1 / (1 + exp(-lambda * (z - .5))), 
  function(z, lambda = c(0.2118881, 0.3006585, 0.0982663, 0.0153671, 0.0016265))
    poly(z, 5) %*% lambda
)

#----- Scenarios

# Number of simulations
ns <- 1000

# Baseline values
N <- 1000
Rho <- 0
Sigma <- .2

# Grid values (other than baseline)
ngrid <- c(50, 100, 200, 500)
rhogrid <- seq(0.25, 0.75, by = 0.25)
sigmagrid <- c(.5, 1)

# Full grid of scenarios
scenarios <- rbind(
  c(n = N, rho = Rho, sigma = Sigma),
  cbind(ngrid, Rho, Sigma),
  cbind(N, rhogrid, Sigma),
  cbind(N, Rho, sigmagrid)
)

#----- Derived objects -----
nsce <- nrow(scenarios)
p <- length(Alpha)
pvec <- sapply(Alpha, length)
ptot <- sum(pvec)

#-------------------------------------------
#     Generate data
#-------------------------------------------

set.seed(1)

# Loop on scenarios
datasim <- apply(scenarios, 1, function(s){
  
  # Correlation matrix
  cormat <- lapply(pvec, function(p){
    sig <- matrix(s["rho"], nrow = p, ncol = p)
    diag(sig) <- 1
    sig
  })
  
  # Generate X
  X <- Map(MASS::mvrnorm, n = s["n"], mu = lapply(pvec, rep, x = 0), 
    Sigma = cormat)
  names(X) <- sprintf("X%i", 1:p)
  
  # Generate indices
  Z <- Map("%*%", X, Alpha)
  
  # Scale Zs
  Z <- lapply(Z, function(z) (z - min(z)) / diff(range(z)))
  
  # Ridge functions
  G <- mapply(function(g, z) do.call(g, list(z = z)), Gfuns, Z)
  
  # Linear predictor
  Y <- 5 + rowSums(scale(G))
  
  # Simulate responses
  Ysim <- Y + matrix(rnorm(s["n"] * ns, 0, s["sigma"]), 
    nrow = s["n"], ncol = ns)
  
  # Return
  c(Y = list(Ysim), X)
})

#-------------------------------------------
#     Apply models
#-------------------------------------------

# Initialize cluster for parallel computation
cl <- makeCluster(max(1, detectCores() - 2))
registerDoParallel(cl)

# To trace simulations
writeLines(c(""), "temp/logsim1.txt")
cat(as.character(as.POSIXct(Sys.time())), file = "temp/logsim1.txt", 
  append = T)

# Packages to load within workers
packs <- c("cgaim", "Matrix", "quadprog", "splines2", "scam")

#----- Start looping on scenarios and simulated responses
results <- foreach(k = seq_len(nsce), .packages = packs, 
    .combine = rbind) %:% 
  foreach(i = seq_len(ns), .packages = packs, .combine = rbind) %dopar%
{
  # Trace iteration
  if(i %% 20 == 0) cat("\n", "scenario =", k, 
    "iter =", i, as.character(Sys.time()), "\n",
      file = "temp/logsim1.txt", append = T)
  
  # Initialize results object
  alpha_est <- list()
  comp_time <- c()
  
  # Current data
  dat <- datasim[[k]]
  dat[[1]] <- dat[[1]][,i]
  
  # Apply GAIM
  comp_time["GAIM"] <- system.time(
    res <- cgaim(Y ~ g(X1) + g(X2) + g(X3), data = dat))[3]
  alpha_est$GAIM <- abs(unlist(res$alpha))

  # Apply CGAIM
  comp_time["CGAIM"] <- system.time(res <- cgaim(Y ~ 
      g(X1, fcons = "inc", acons = list(monotone = -1, sign = 1)) + 
      g(X2, fcons = "inc", acons = list(monotone = 1, sign = 1)) + 
      g(X3, fcons = "cvx", acons = list(sign = 1)),
    data = dat, smooth_control = list(sp = rep(0, p))))[3]
  alpha_est$CGAIM <- unlist(res$alpha)

  # Apply misspecified GAIM
  comp_time["MGAIM"] <- system.time(res <- cgaim(Y ~ 
      g(X1, fcons = "inc", acons = list(monotone = 1, sign = 1)) + 
      g(X2, fcons = "inc", acons = list(monotone = -1, sign = 1)) + 
      g(X3, fcons = "ccv", acons = list(sign.const = 1)),
    data = dat, smooth_control = list(sp = rep(0, p))))[3]
  alpha_est$MGAIM <- unlist(res$alpha)

  # Apply cumulative effect GAM (Xia & Tong, 2005)
  x <- datasim[[k]][-1]; x[[1]] <- x[[1]][,pvec[1]:1]
  comp_time["FACTS"] <- system.time(
    res <- facts(x, dat$Y, gshape = c("inc", "inc", "cvx"), 
      thetashape = c("inc", "inc", "nc")))[3]
  alpha <- res$theta; alpha[[1]] <- rev(alpha[[1]])
  alpha_est$FACTS <- unlist(alpha)

  # Apply gMAVE
  comp_time["gMAVE"] <- system.time(
    res <- gmave(dat$Y, datasim[[k]][-1], alpha.norm = "1")
  )[3]
  alpha_est$gMAVE <- abs(res)

  # Apply PPR
  comp_time["PPR"] <- system.time(
    res <- ppr_apply(datasim[[k]][-1], dat$Y, norm.type = "1")
  )[3]
  alpha_est$PPR <- abs(res)
  
  # Return
  data.frame(scenarios[k,,drop = F], simu = i, 
    alpha = rep(1:ptot, length(alpha_est)),
    model = rep(names(alpha_est), each = ptot), 
    est = unlist(alpha_est, use.names = F),
    time = rep(comp_time, each = ptot))
}

# Close workers
stopCluster(cl)

#-------------------------------------------
#  Prepare plots
#-------------------------------------------

# Palette and symbols
pal <- c(rev(brewer.pal(5, "Blues")[-(1:2)]), 
  colorRampPalette(c("saddlebrown", "white"))(5)[-(4:5)])
pchs <- rep(16:17, each = 3)
names(pal) <- names(pchs) <- 
  c("CGAIM", "GAIM", "MGAIM", "gMAVE", "FACTS", "PPR")

#-------------------------------------------
#  Compute RMSE
#-------------------------------------------

# Get number of models
nm <- length(unique(results$model))

# Add true Alpha
results$true <- rep(unlist(Alpha), nsce * ns * nm)

# Compute Squared-error
results$sqerr <- with(results, (est - true)^2)

# Compute RMSE of specific simulation
rmse <- aggregate(sqerr ~ n + rho + sigma + model + simu, data = results,
  function(x) sqrt(mean(x)))

# Compute the average and confidence intervals
err_summary <- aggregate(sqerr ~ n + rho + sigma + model, data = rmse, 
  function(x) c(mean(x), sd(x) / ns))
colnames(err_summary[["sqerr"]]) <- c("est", "se")
err_summary <- cbind(err_summary[names(err_summary) != "sqerr"], 
  err_summary[["sqerr"]])
err_summary$high <- with(err_summary, est + 1.96 * se)
err_summary$low <- with(err_summary, est - 1.96 * se)

#----- Figure 1: RMSE

# Create common plot layout
plotlay <- ggplot(err_summary, aes(y = est, ymin = low, 
    ymax = high, col = model, fill = model, shape = model)) + 
  theme_classic() + 
  geom_ribbon(alpha = .2, col = NA) +
  geom_line(size = 1) + 
  geom_point(size = 3) +
  scale_fill_manual(values = pal, name = "Method") + 
  scale_color_manual(values = pal, name = "Method") + 
  scale_shape_manual(values = pchs, name = "Method") +
  scale_y_continuous(name = "RMSE", trans = "log10") + 
  coord_cartesian(ylim = c(4e-3, .4)) +
  theme(panel.grid.major.y = element_line(colour = "grey", linetype = 2))

# Plot by n
plotn <- plotlay %+% subset(err_summary, rho == Rho & sigma == Sigma) + 
  aes(x = n) + xlab(expression(n)) + ggtitle("a) Sample size")

# Plot by sigma
plotsigma <- plotlay %+% subset(err_summary, rho == Rho & n == N) + 
  aes(x = sigma) + xlab(expression(sigma)) + ggtitle("c) Noise level")

# Plot by rho
plotrho <- plotlay %+% subset(err_summary, sigma == Sigma & n == N) + 
  aes(x = rho) + xlab(expression(rho)) + ggtitle("b) Correlation")

# Put together
design <- "
  1122
  #33#
"
plotn + plotrho + plotsigma + 
  plot_layout(design = design, guides = "collect", 
    widths = rep(1, 4), heights = rep(1,2))

#-------------------------------------------
#  Detail of alphas
#-------------------------------------------

# Transform model as factor
results$model <- factor(results$model, 
  levels = c("CGAIM", "GAIM", "MGAIM", "gMAVE", "FACTS", "PPR"))

# Supplementary Figure 3
alphaplot <- ggplot(results) + theme_classic() +
  geom_boxplot(aes(x = factor(alpha), y = est, color = model, shape = model)) + 
  scale_color_manual(values = pal, name = "Method") + 
  scale_shape_manual(values = pchs, name = "Method") + 
  geom_segment(aes(x = alpha - .4, xend = alpha + .4, y = true, yend = true)) +
  geom_vline(xintercept = 1:(ptot -1) + .5, col = "grey", linetype = 2) +
  geom_vline(xintercept = cumsum(pvec[-length(pvec)]) + .5) +
  scale_x_discrete(name = "Variable", 
    labels = c(sapply(pvec, seq, from = 1, by = 1))) + 
  scale_y_continuous(name = "Estimated", limits = c(0, 1), expand = c(0, 0)) + 
  annotate("label", x = tapply(1:ptot, rep(1:p, pvec), mean), y = 0.95, 
    label = sprintf("Index %i", 1:p), size = 5, 
    fill = "white") + 
  theme(plot.margin = unit(c(3, 1, 1, 1), "lines"), 
    plot.title = element_text(hjust = .5, size = 15),
    panel.border = element_rect(fill = NA, color = "black"))

# Execute plot for each scenario
allplots <- lapply(seq_len(nsce), function(j) alphaplot %+% subset(results,
  n == scenarios[j,"n"] & rho == scenarios[j,"rho"] & 
    sigma == scenarios[j,"sigma"]) + 
  ggtitle(bquote(n == .(scenarios[j,"n"]) * "," ~ 
      rho == .(scenarios[j,"rho"]) * "," ~ sigma == .(scenarios[j,"sigma"]))))

# Put together
wrap_plots(allplots, ncol = 2, nrow = 5, guides = "collect")

#-------------------------------------------
# Computation time
#-------------------------------------------

# Extract summary
time_sum <- aggregate(time ~ model, data = results, summary)

# Put Mean and IQR in a string
char_time <- apply(time_sum[[2]], 1, function(x) 
  sprintf("%2.2f (%2.2f - %2.2f)", x["Mean"], x["1st Qu."], x["3rd Qu."]))

# Supplementary Table 4
time_tab <- cbind(Model = time_sum[[1]], Time = char_time)
time_tab <- time_tab[match(names(pal), time_tab[,1]),]
