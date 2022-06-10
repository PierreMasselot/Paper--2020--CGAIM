################################################################################
#
#   Wrapper for the PPR used in simulation in
#
#   Masselot et al., 2022
#   Constrained groupwise additive index models
#   Biostatistics
#
#   Author: Pierre Masselot
#
################################################################################

#----- Parameters
# x: list of group matrices
# y: response vector
# norm.type: norm constraining parameters. Same as in cgaim package.
# ...: parameters to be passed to ppr function

#----- Value
# A vector of stacked weights alpha. Follows the order of variables in x.

ppr_apply <- function(x, y, norm.type = "L2",  ...)
{
  # Prepare data
  x <- lapply(x, as.matrix)
  p <- length(x)
  pvec <- sapply(x, ncol)
  d <- sum(pvec)
  pind <- rep(1:p, pvec)
  Xall <- Reduce(cbind, x)
  nx <- sapply(x, nrow)
  ny <- length(y)
  n <- unique(c(nx, ny))
  if (length(n) != 1){
    stop("inconsistent number of observation")
  }
  # Apply ppr
  res <- ppr(y = y, x = Xall, nterms = p, ...)
  # Extract alphas
  alpha <- Map("[", as.data.frame(res$alpha), 
    split(seq_len(ptot), rep(1:p, pvec)))  # Alphas
  unlist(Map(cgaim:::normalize, alpha, norm.type))
}