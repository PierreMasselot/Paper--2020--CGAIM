###############################################################################
#
#                          Simulation study
#                          Plotting results
#
###############################################################################

#-------------------------------------------
#    Parameters
#-------------------------------------------

true_alphas <- c(.7, .2, .1, 0, 0, 0, .5, .5, .2, .4, .3, .1)

# Graphical parameters
mod_names <- c("GAIM", "CGAIM", "PPR", "gMAVE")
mod_col <- c(4, 5, 2, 3)
mod_pch <- 15:18
mod_lty <- 

par_names <- c("Sample size", "Correlation coefficient", "Sample size", 
  "Sample size")
exp_names <- c("CGAIM structure", "Correlated variables", 
  "Misspecified model", "No group structure")

#-------------------------------------------
#    Loading data
#-------------------------------------------

# List all result files
fillist <- list.files("Results", pattern = "1[.]")
nfil <- length(fillist)

# Prepare result storing objects
alpha_errors <- vector("list", nfil)
par_vecs <- vector("list", nfil)
for (i in seq_along(fillist)){
  # Load result
  objs <- load(sprintf("Results/%s", fillist[i]))
  # Compute errors of estimated alphas
  aerr <- lapply(alpha_samp, lapply, apply, 2, "-", true_alphas)
  # Compute mean of errors
  alpha_errors[[i]] <- sapply(aerr, sapply, apply, 1, 
    function(x) sqrt(sum(x^2, na.rm = T)), simplify = "array")
  # Extract the parameter vector object
  parobj <- grep("vec", objs)
  par_vecs[[i]] <- get(objs[parobj])
}

#-------------------------------------------
#  Plots  
#-------------------------------------------
pl_dim <- n2mfrow(nfil)
lay_mat <- rbind(matrix(seq_len(nfil), pl_dim[1], pl_dim[2], byrow = T),
  nfil + 1)

x11(width = 20, height = 15)
layout(lay_mat, heights = c(rep(1, pl_dim[1]), .2))
for (i in seq_len(nfil)){
  alph_errs <- apply(alpha_errors[[i]], 2:3, mean)
  
  matplot(par_vecs[[i]], alph_errs, type = "b", log = "y", col = mod_col,
    pch = mod_pch, xlab = par_names[[i]], ylab = "RMSE", lwd = 2, cex = 1.5,
    ylim = c(.1, 10000), cex.lab = 1.5, cex.axis = 1.3,
    main = sprintf("%s) %s", letters[i], exp_names[i]), cex.main = 1.5)
}
par(mar = rep(0,4))
plot.new()
legend("center", mod_names, col = mod_col, lwd = 2, lty = 1:5, 
  ncol = length(mod_names), cex = 1.5, pch = mod_pch)

dev.print(png, filename = "Results/Figure1.png", res = 200, 
  width = dev.size()[1], height = dev.size()[2], units = "in")
