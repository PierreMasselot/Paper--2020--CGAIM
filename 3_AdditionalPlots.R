################################################################################
#
#  Additional graphics found in the Supplementary Materials of
#
#   Masselot et al., 2022
#   Constrained groupwise additive index models
#   Biostatistics
#
#   Supplementary materials
#
#   Author: Pierre Masselot
#
################################################################################

#-------------------------------------------
# Packages
#-------------------------------------------

library(viridis) # Color palettes
library(corrplot) # To plot correlation matrix

#-----------------------
# Supplementary Figure 1
#-----------------------

# True functions
Gfuns <- c(
  function(z) log(z - min(z) + .1), 
  function(z, lambda = 20) 1 / (1 + exp(-lambda * (z - .5))), 
  function(z, lambda = c(0.2118881, 0.3006585, 0.0982663, 0.0153671, 0.0016265))
    poly(z, 5) %*% lambda
)

# Evaluate functions
zseq <- seq(0, 1, length.out = 1000)
geval <- sapply(Gfuns, do.call, list(z = zseq))

# Plot
matplot(zseq, scale(geval), type = "l", lty = 1, lwd = 2, 
  col = viridis(length(Gfuns)), xlab = "Standardized index", ylab = "g(.)")
legend("topleft", legend = sapply(seq_along(Gfuns), 
    function(i) as.expression(bquote(g[.(i)](.)))),
  col = viridis(length(Gfuns)), lty = 1, lwd = 2, bty = "n")

#-----------------------
# Supplementary Figure 2
#-----------------------

# Load correlation matrix (Robinson et al. 2018)
corrmat <- data.matrix(read.table("Data/1.4_ExposomeCorrelation.csv", 
  sep = ","))
colnames(corrmat) <- rownames(corrmat) <- NULL

# Group specification
groups <- list(meteorology = c(1:4), ap = c(5:9), road = c(10, 16:18), 
  natural = c(11:15), built = c(19:28))
groupinds <- unlist(groups)
groupvec <- rep(seq_along(groups), sapply(groups, length))

# Reorder correlation matrix
corrmat_ord <- corrmat[groupinds, groupinds]

# Limits of groups
grouplims <- which(diff(groupvec) > 0) + .5

# Plot correlation matrix
corrplot(corrmat_ord, tl.pos = "n")
segments(x0 = grouplims, y0 = .5, y1 = length(groupinds) + .5)
segments(x0 = .5, x1 = length(groupinds) + .5, 
  y0 = length(groupinds) - grouplims + 1)
rect(.5, length(groupinds) + .5, length(groupinds) + .5, .5, 
  col = NA, border = "black")
