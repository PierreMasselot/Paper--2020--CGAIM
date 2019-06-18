######################################################
#
#              Simulation results
#
######################################################
library(RColorBrewer)

source("C:/Users/masselpl/Documents/Recherche/# Ressources Perso/R/Fonctions/Plot_functions.R")
source("C:/Users/masselpl/Documents/Recherche/# Ressources Perso/R/Fonctions/Util_functions.R")

setwd("C:/Users/masselpl/Documents/Recherche/2017-2019 - Post-doc/Programmes/R/2 - Indices PPR/Paper--2020--CGAIM")

#----------------------------------
#          Parameters
#----------------------------------

models <- c("GNConst", "PPR", "MAVE")
mod_names <- c("cGAIM", "PPR", "RMAVE")

nm <- length(models)

#----------------------------------
#       Loading results
#----------------------------------

load("Simulation_results.RData")

final_mod <- mod[models]

#----------------------------------
#         Alphas
#----------------------------------

# Show Alphas

x11(title = "Alpha_coefficients")
par(mfrow = c(1,p), mar = c(5.1,0,0,0), oma = c(0, 5.1, 4.1, 1.1))
for (j in 1:p){
  bplist <- vector("list", nm * pvec[j])
  cp <- 0
  for (l in 1:pvec[j]){
    for (i in 1:nm){
      cp <- cp + 1
      bplist[[cp]] <- sapply(final_mod[[i]], function(x) x$alpha[[j]][l])
      bplist[[cp]] <- bplist[[cp]] - Alpha[[j]][l]
    }
  }
  at <- 1:cp + 0:(cp - 1) %/% nm
  cuts <- setdiff(1:max(at), at)
  boxplot(bplist, at = at, xaxt = "n", border = brewer.pal(nm,"Set1"), lwd = 2,
    xlab = substitute(expression(alpha[j]), list(j = j)), cex.axis = 1.3, 
    cex.lab = 1.5, ylim = c(-1.5,1.5), yaxt = "n")
  if (j == 1){
    axis(2, xpd = T, cex.axis = 1.2)
    title(ylab = "Alpha error", cex.lab = 1.5, xpd = NA)
    outerLegend("topleft", mod_names, lwd = 2, ncol = 3, cex = 1.5,
      col = brewer.pal(nm,"Set1"), bty = "n")
  } 
  abline(v = cuts, lty = 2, col = "grey")
  abline(h = 0, lwd = 3)
  axis.intervals(ticks = c(par("usr")[1], cuts, par("usr")[2]), 
    labels = 1:pvec[j] - 1, cex.axis = 1.3)
}

dev.print(png, filename = "Alpha_coeffcients.png", res = 100, 
  width = dev.size()[1], height = dev.size()[2], units = "in")
  

# Show errors
alphas_diff <- lapply(final_mod, lapply, function(x){
  mapply("-", x$alpha, Alpha)
})
weight_SE <- sapply(alphas_diff, sapply, function(x)
  unlist(x)^2
)
overall_RMSE <- sqrt(apply(weight_SE, 2, mean))

x11("overall_RMSE")
barplot(overall_RMSE, col = brewer.pal(nm,"Set1"), border = NA, log = "y",
  ylab = "RMSE (log)", cex.lab = 1.3, cex.axis = 1.2, cex.names = 1.2, 
  names.arg = mod_names)

dev.print(png, filename = "overall_RMSE.png", res = 100, 
  width = dev.size()[1], height = dev.size()[2], units = "in")

  
# Ridge functions
x11(title = "True_functions", width = 10, height = 4)
par(mfrow = c(1, p), mar = c(2, 2, 4, 2))
for (j in 1:p){
  zord <- order(Z[,j])
  plot(Z[zord,j], G[zord,j], type = "l", lwd = 3, xaxt = "n", yaxt = "n", 
    main = substitute(g[j](), list(j = j)), cex.main = 1.5)
}

dev.print(png, filename = "True functions.png", res = 100, 
  width = dev.size()[1], height = dev.size()[2], units = "in")


# Estimated ridge functions
zseq <- seq(0, 1, length.out = 1000)
gs <- vector("list", p)
for (j in 1:p){
  gint <- lapply(final_mod, sapply, function(x){
    spline(x$z[,j], x$gz[,j], xout = zseq)$y
  })
  gs[[j]] <- sapply(gint, apply, 1, mean)
}

  
x11(title = "estimated_functions", width = 10, height = 4)
par(mfrow = c(1, p))
for (j in 1:p){
  gs[[j]][,3] <- predict(loess(gs[[j]][,3]~zseq))
  zord <- order(Z[,j])
  matplot(zseq, scale(gs[[j]]), type = "l", lwd = 2, lty = 1, 
    col = brewer.pal(nm,"Set1"), xlab = substitute(I[j](), list(j = j)), 
    ylab = substitute(g[j](), list(j = j)), cex.lab = 1.3, cex.axis = 1.2)
  lines(fscaling(Z[zord,j]), scale(G[zord,j]), lwd = 3)
  if (j == 1) outerLegend("topleft", mod_names, lwd = 2, ncol = 3, cex = 1.5,
    col = brewer.pal(nm,"Set1"), bty = "n")
}


dev.print(png, filename = "Estimated_functions.png", res = 100, 
  width = dev.size()[1], height = dev.size()[2], units = "in")