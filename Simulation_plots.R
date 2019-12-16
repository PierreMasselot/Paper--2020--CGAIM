######################################################
#
#              Simulation results
#
######################################################
library(RColorBrewer)
library(sfsmisc)

source("C:/Users/masselpl/Documents/Recherche/# Ressources Perso/R/Fonctions/Plot_functions.R")
source("C:/Users/masselpl/Documents/Recherche/# Ressources Perso/R/Fonctions/Util_functions.R")

setwd("C:/Users/masselpl/Documents/Recherche/2017-2019 - Post-doc/Resultats/Part 2 - GAIM/Article V3")

#----------------------------------
#          Parameters
#----------------------------------

mod_names <- c("GAIM", "CGAIM", "MGAIM", "PPR", "gMAVE")
nm <- length(mod_names)

mod_pal <- c("royalblue3", "cornflowerblue", "lightblue1", "indianred", 
  "forestgreen")
mod_pch <- 14 + 1:nm
mod_lty <- 1:nm

# Single design
Exp <- 4

varPar_name = "Correlation coefficient"

SAVE <- FALSE

#----------------------------------
#       Loading results
#----------------------------------

load("C:/Users/masselpl/Documents/Recherche/2017-2019 - Post-doc/Resultats/Part 2 - GAIM/Article V3/Collinearity_study_resultsV3.RData")

nexp <- length(g_all[[1]])
p <- length(constant_parameters$Alpha)
pvec <- sapply(constant_parameters$Alpha, length)
pind <- rep(1:p, sapply(constant_parameters$Alpha, length))
ns <- constant_parameters$ns

#----------------------------------
#         Alphas
#----------------------------------

true_alphas <- unlist(constant_parameters$Alpha)
  
# Alphas of a single design
exp_alphas <- lapply(alpha_all, "[[", Exp)

x11(title = "Alpha_SingleDesign", height = 5, width = 10)
layout(rbind(1:p, p + 1), heights = c(.9, .1))
for (j in 1:p){
  bpdat <- list()
  for (i in which(pind == j)){
    bpdat <- c(bpdat, lapply(exp_alphas, "[", i, ))
  } 
  at <- 1:(pvec[j]*nm) + 0:((pvec[j]*nm) - 1) %/% nm
  cuts <- setdiff(1:max(at), at)
  bplims <- sapply(bpdat, function(x) boxplot.stats(x)$stats[c(1,5)])
  boxplot(bpdat, at = at, xaxt = "n", border = mod_pal, lwd = 2,
    ylab = substitute(expression(alpha[j]), list(j = j)), cex.axis = 1.3, 
    cex.lab = 1.5, ylim = range(c(unlist(bplims), -1, 1)), outpch = ".")
  abline(v = cuts, lty = 2, col = "grey")
  segments(c(0, cuts), true_alphas[pind == j], 
    c(cuts, max(at) + 1), true_alphas[pind == j], lwd = 3)
  axis.intervals(ticks = c(par("usr")[1], cuts, par("usr")[2]), 
    labels = 1:pvec[j], cex.axis = 1.3)
}
par(mar = rep(0,4))
plot.new()
legend("center", mod_names, fill = mod_pal, cex = 2, border = NA,
  ncol = nm, bty = "n")
if (SAVE){  
  dev.print(png, filename = "Simu_alphasOneCase.png", res = 100, 
    width = dev.size()[1], height = dev.size()[2], units = "in")
}
  

#----------------------------------
#          Final plots
#----------------------------------

# Functions error (MISE)

trueGs <- lapply(datSim, "[[", "G")
trueZs <- lapply(datSim, function(x){
  mapply("%*%", x$X, constant_parameters$Alpha)
})

ISEs <- array(NA, dim = c(nexp, p, nm, ns))
for (e in 1:nexp){
  for (j in 1:p){
    trueFunc <- splinefun(trueZs[[e]][,j], trueGs[[e]][,j])
    trueEval <- trueFunc(seq(min(trueZs[[e]][,j]), max(trueZs[[e]][,j]), 
      length.out = 1000))
    for (m in 1:nm){
      modFunc <- Map(function(z, g){
          if (all(is.na(g))){
            out <- rep(NA, 1000)
          } else{
            out <- splinefun(z, g)(seq(min(z), max(z), length.out = 1000))
          }
          return(out)
        }, 
        as.data.frame(z_all[[m]][[e]][,j,]),
        as.data.frame(g_all[[m]][[e]][,j,]))
      # Integrate the difference by trapezoid approximation
      ISEs[e,j,m,] <- mapply(function(z, g){
          ifelse(all(is.na(g)) || all(g == 0), NA, 
            integrate.xy(seq(min(z), max(z), length.out = 1000), 
              (g - trueEval)^2))
        }, as.data.frame(z_all[[m]][[e]][,j,]), modFunc) 
    }
  }
}
MISEs <- apply(ISEs, 1:3, mean, na.rm = T)


# Alpha errors
true_alphas <- unlist(constant_parameters$Alpha)

# Compute RMSEs for each alpha
alpha_errors <- lapply(alpha_all, lapply, apply, 2, "-", true_alphas)
alpha_mses <- lapply(alpha_errors, sapply, apply, 1, 
  function(x) sqrt(sum(x^2, na.rm = T)))

alpha_mse_ses <- Map(function(er, mse){
  mapply(function(er1, mse1){
    vars <- apply(er1, 2, function(x) x^2 - mse1)
    rowSums(vars) / (ns * (ns - 1))
  }, er, as.data.frame(mse))
}, alpha_errors, alpha_mses)

alpha_mses_byZ <- lapply(alpha_mses, aggregate, by = list(Group = pind), mean)
error_rg <- range(sapply(alpha_mses, function(x) range(x[,-1])))

# plot
x11(title = "Simu_results", height = 10, width = 15)
layout(rbind(matrix(1:(2*p), nrow = 2), (2*p) + 1), 
  heights = c(.45, .45, .1))
for (j in 1:p){
  group_errors <- sapply(alpha_mses_byZ, "[", j, -1)
  xpars <- unlist(varying_parameters)
  matplot(xpars, group_errors, type = "b", lwd = 3,
    col = mod_pal, pch = mod_pch, log = "y", xlab = varPar_name,
    main = bquote(.(letters[j]) * ") Group" ~ .(j) * ":" ~ italic(alpha)),
      #paste0(letters[j], ") Group ", j, expression(alpha)),
    ylim = error_rg, 
    ylab = ifelse(j == 1, "MSE", ""), cex.lab = 1.5, cex.axis = 1.3, 
    cex.main = 2, cex = 1.5, lty = mod_lty)
  matplot(xpars, MISEs[,j,], type = "b", lwd = 3,
    col = mod_pal, pch = mod_pch, log = "y", xlab = varPar_name,
    main = bquote(.(letters[p + j]) * ") Group" ~ .(j) * ":" ~ italic(g)), 
    ylim = range(MISEs), 
    ylab = ifelse(j == 1, "MISE", ""), cex.lab = 1.5, cex.axis = 1.3, 
    cex.main = 2, cex = 1.5, lty = mod_lty)
}
par(mar = rep(0,4))
plot.new()
legend("center", mod_names, col = mod_pal, lwd = 2, lty = mod_lty, 
  ncol = nm, cex = 1.5, pch = mod_pch)

dev.print(png, filename = sprintf("Simu_result_%s.png", varPar_name), res = 100, 
  width = dev.size()[1], height = dev.size()[2], units = "in")