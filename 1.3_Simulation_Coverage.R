################################################################################
#
#  R code for the simulation study of 
#
#   Masselot et al., 2022
#   Constrained groupwise additive index models
#   Biostatistics
#
#   Section 4.3 Coverage
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
library(data.table) # For function between

#----- cgaim package
# Should be installed from github
# install_github("PierreMasselot/cgaim")
library(cgaim)

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
  function(z, lambda = c(0.2118881, 0.3006585, 0.0982663, 0.0153671, 0.0016265))
    poly(z, 5) %*% lambda
)

# Sample size
n <- 1000

# Number of simulations
ns <- 1000
B <- 500

# Noise levels
Sigmas <- c(.2, .5, 1)

#----- Derived objects -----
p <- length(Alpha)
pvec <- sapply(Alpha, length)
ptot <- sum(pvec)

#-------------------------------------------
#     Simulations
#-------------------------------------------

#----- Generate data    
set.seed(3)

X <- Map(MASS::mvrnorm, n = n, mu = lapply(pvec, rep, x = 0), 
  Sigma = lapply(pvec, diag))
names(X) <- sprintf("X%i", 1:p)
Z <- Map("%*%", X, Alpha)
G <- mapply(function(g, z) do.call(g, list(z = z)), Gfuns, Z)
Y <- 5 + rowSums(scale(G))
Ysim <- lapply(Sigmas, function(s) Y + matrix(rnorm(n * ns, 0, s), n, ns))

#----- Alpha estimation by all algorithms

# Initialize cluster for parallel computation
cl <- makeCluster(max(1, detectCores() - 2))
registerDoParallel(cl)

# To trace simulations
writeLines(c(""), "temp/logsim3.txt")
cat(as.character(as.POSIXct(Sys.time())), file = "temp/logsim3.txt", 
  append = T)

# Packages
packs <- c("cgaim", "Matrix", "scar", "quadprog", "splines2", "scam", "foreach")

# foreach loops
results <- foreach(s = seq_along(Sigmas), .packages = packs, 
    .combine = rbind) %:% 
  foreach(i = seq_len(ns), .packages = packs, 
    .combine = rbind) %dopar% 
{
  # Trace iteration
  if(i %% 20 == 0) cat("\n", "scenario =", s, 
    "iter =", i, as.character(Sys.time()), "\n",
      file = "temp/logsim3.txt", append = T)
  
  # Initialize data and result object
  y <- Ysim[[s]][,i]
  dat <- c(list(y = y), X)

  # Apply CGAIM
  res <- cgaim(y ~ 
      g(X1, fcons = "inc", acons = list(monotone = -1, sign = 1)) + 
      g(X2, fcons = "inc", acons = list(monotone = 1, sign = 1)) + 
      g(X3, fcons = "cvx", acons = list(sign = 1)),
    data = dat, smooth_control = list(sp = rep(0, p)))
  
  #----- Compute confidence intervals
  cis <- list()
  
  # Compute Bootstrap confidence interval by residual resampling
  bootres <- boot.cgaim(res, B = B, boot.type = "residuals")
  cis$Bootstrap <- confint(bootres, parm = 1)$alpha
  
  # Compute normal confidence interval
  cis$Normal <- confint(res, B = B, type = "norm", parm = 1)$alpha
  
  # Put them all together
  allcis <- do.call(rbind, cis)
  colnames(allcis) <- c("low", "high")
  
  # output
  data.frame(sigma = Sigmas[s], simu = i, 
    type = rep(names(cis), each = ptot),
    alpha = rep(1:ptot, length(cis)), 
    est = rep(unlist(res$alpha), length(cis)), 
    allcis)
}

# Stop parallel
stopCluster(cl)

#-------------------------------------------
#     Results
#-------------------------------------------

#----- Estimate overall coverage

# Check if includes true alpha
eps <- alpha.control()$ctol
results$includes <- between(unlist(Alpha)[results$alpha], results$low - eps,
  results$high + eps)

# Compute average estimated alpha
mean_est <- aggregate(est ~ alpha + type + sigma, results, mean)
names(mean_est)[length(mean_est)] <- "mean_est"
results <- merge(results, mean_est)

# Check if include mean estimated alpha
results$includes_mean <- between(results$mean_est, results$low, results$high)

# Compute coverage for each
ov_coverage <- aggregate(cbind(includes, includes_mean) ~ type + sigma, 
  results, function(x) c(coverage = mean(x), 
    se = sqrt(mean(x) * (1 - mean(x)) / ns)))
ov_coverage <- do.call(data.frame, ov_coverage)
names(ov_coverage)[-(1:2)] <- apply(
  expand.grid(c("est", "se"), c("cover", "bc_cover")), 1, 
  function(x) paste(rev(x), collapse = "_"))

#----- Estimate alpha-specific coverage

# Compute coverage
alpha_coverage <- aggregate(
  cbind(includes, includes_mean) ~ type + alpha + sigma, data = results, 
  function(x) c(coverage = mean(x), se = sqrt(mean(x) * (1 - mean(x)) / ns)))
alpha_coverage <- do.call(data.frame, alpha_coverage)
names(alpha_coverage)[-(1:3)] <- apply(
  expand.grid(c("est", "se"), c("cover", "bc_cover")), 1, 
  function(x) paste(rev(x), collapse = "_"))

#-------------------------------------------
#     Plots
#-------------------------------------------

# Number of confidence interval types
ntypes <- length(unique(ov_coverage$type))

# Color palette
colpal <- tail(brewer.pal(pmax(ntypes, 3), "Blues"), ntypes)

#----- Figure 3: Overall bias-corrected coverage

ggplot(ov_coverage) + theme_classic() +
  geom_pointrange(aes(x = factor(sigma), y = bc_cover_est, 
    ymin = bc_cover_est - bc_cover_se, 
    ymax = bc_cover_est + bc_cover_se, 
    group = type, col = type, shape = type), 
    position = position_dodge(width = .2)) + 
  geom_hline(yintercept = c(.95), linetype = 2) + 
  labs(y = "Coverage", shape = "Type", x = "Noise level") + 
  scale_color_manual(name = "Type", values = colpal) +
  theme(axis.ticks.x.top = element_line(linetype = 0),
    axis.text.x.top = element_text(size = 12),
    axis.line.x.top = element_blank(),
    panel.border = element_rect(fill = NA),
    panel.grid.major = element_line(color = "grey", linetype = 3)) + 
  ylim(c(min(with(ov_coverage, bc_cover_est - bc_cover_se)), 1))

#----- Supplementary Figure 4: Alpha specific bias-corrected coverage 

ggplot(alpha_coverage) + theme_classic() +
  geom_pointrange(aes(x = alpha, y = bc_cover_est, 
    ymin = bc_cover_est - bc_cover_se, 
    ymax = bc_cover_est + bc_cover_se, 
      group = type, col = type, shape = type), 
    position = position_dodge(width = .5)) + 
  facet_wrap(~ sigma, labeller = label_bquote(sigma == .(sigma))) +
  geom_hline(yintercept = c(.95), linetype = 2) + 
  geom_vline(xintercept = cumsum(pvec) + .5, linetype = 2, col = "grey") + 
  labs(y = "Coverage", shape = "Type") + 
  scale_color_manual(name = "Type", 
    values = brewer.pal(ntypes + 1, "Blues")[-1]) +
  scale_x_continuous(name = "Covariate", breaks = seq_len(ptot), 
    labels = unlist(lapply(pvec, seq, from = 1)),
    sec.axis = sec_axis(trans = ~., name = "",
      breaks = tapply(1:ptot, rep(1:p, pvec), mean),
      labels = sprintf("Index %i", 1:p))) + 
  theme(axis.ticks.x.top = element_line(linetype = 0),
    axis.text.x.top = element_text(size = 12),
    axis.line.x.top = element_blank(),
    panel.border = element_rect(fill = NA),
    panel.grid.major.y = element_line(color = "grey", linetype = 3))
