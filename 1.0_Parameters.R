###############################################################################
#
#                          Simulation study
#                             Parameters
#
###############################################################################


#---- Simulation parameters
ns <- 5000 # Number of simulations
n <- 1000

Beta0 <- 0  # Intercept of the whole model
Beta1 <- rep(1 / 3, 3) # Index importance
Alpha <- list(c(.7, .2, .1), c(.2, .3, .5), c(.2, .6, .2)) # Index weights

Gfuns <- c("g1", "g2", "g3")  # Ridge functions
Gpars <- list(list(lambda = 1), 
  list(lambda = 5), 
  list(a1 = 50, a2 = 0.8, b1 = .5, b2 = 2.3)
) # Parameters for Ridge functions

Xcorr <- 0

Ysigma <- .2    # Noise amplitude

#---- Plot parameters

mod_names <- c("GAIM", "CGAIM", "MGAIM", "PPR", "gMAVE")
mod_pal <- c("royalblue3", "cornflowerblue", "lightblue1", "indianred", 
  "forestgreen")
mod_pch <- 15:19
mod_lty <- 1:5

#---- useful quantities
p <- length(Alpha)
pvec <- sapply(Alpha, length)
d <- sum(pvec)
pind <- rep(1:p, pvec)
nmod <- length(mod_names)