library(splines)
library(mgcv)
library(devtools)
library(scam)
library(forecast)
library(parallel)

setwd("C:/Users/masselpl/Documents/Recherche/2017-2019 - Post-doc/Programmes/R/2 - Indices PPR/Paper--2020--CGAIM")

# Should be installed from github
# install_github("PierreMasselot/cgaim")
library(cgaim)
source("Useful_functions.R")

respath <- "C:/Users/masselpl/Documents/Recherche/2017-2019 - Post-doc/Resultats/Part 2 - GAIM/Article V4"

#---------------------------------------------------
#  Data reading
#---------------------------------------------------

datatab <- read.table("Data2 - Pollution.csv", sep = ";", header = T)

# Compute temperature 3-day mean
datatab$T3d <- ma(datatab$Tmean, 3)
datatab$numday <- as.numeric(datatab$Date)

datatab <- na.omit(datatab)
n <- nrow(datatab)

#---------------------------------------------------
#     Applying CGAIM
#---------------------------------------------------

cgaim_res <- cgaim(MCV ~ 
  g(NO2, O3, PM2.5, bs = "mpi", constraints = list(sign.const = 1)) + 
  s(numday) + s(T3d), data = datatab,
  alpha.control = list(norm.type = "sum"),
  algo.control = list(keep.trace = T))

p <- ncol(cgaim_res$gfit)
d <- ncol(cgaim_res$indexfit)

#---------------------------------------------------
#    Bootstrapping for uncertainty evaluation
#---------------------------------------------------

cl <- makeCluster(2)
cgaim_CI <- confint(cgaim_res, B = 4, l = 7, 
  applyFun = "parLapply", cl = cl)
stopCluster(cl)

if (F) save(cgaim_res, cgaim_CI, 
  file = sprintf("%s/Bootstrap_results_2.RData", respath))
if (F) load(sprintf("%s/Bootstrap_results_2.RData", respath))
  
#---------------------------------------------------
#     Visualizing results
#---------------------------------------------------

bca_cgaim <- cgaim_CI$alpha$boot.bca

x11(width = 10, height = 5)
par(mfrow = c(1,2), cex.lab = 1.3, cex.axis = 1.2)
  
plot(cgaim_res$alpha[[1]], ylim = c(0,1), col = "white",
  xlab = "Pollutant", ylab = expression(alpha), xaxt = "n",
  main = bquote("a) Index weights" ~ italic(alpha)))
axis(1, at = 1:p, 
  labels = c(expression(NO[2]), expression(O[3]), expression(PM[2.5])))
arrows(x0 = 1:p, y0 = bca_cgaim[,1], y1 = bca_cgaim[,2],
  angle = 90, length = 0.05, lwd = 2, code = 3, col = "darkgrey")
points(1:p, cgaim_res$alpha[[1]], col = "darkgrey", pch = 21, cex = 2, 
  bg = "cornflowerblue")
  
plot(cgaim_res, select = 1, xlab = "Z", ylab = "g()", 
  col = "cornflowerblue", lwd = 3,
  main = bquote("b) Function" ~ italic(g())),
  ci = cgaim_CI, ci.args = list(col = "lightgrey")
)

dev.print(png, filename = sprintf("%s/App2_Results.png", respath), res = 60, 
  width = dev.size()[1], height = dev.size()[2], units = "in")
dev.copy2eps(file = sprintf("%s/App2_Results.eps", respath))
