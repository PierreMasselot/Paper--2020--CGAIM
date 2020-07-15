library(splines)
library(mgcv)
library(devtools)
library(scam)
library(parallel)
library(forecast)


# Should be installed from github
# install_github("PierreMasselot/cgaim")
library(cgaim)
source("0_Useful_functions.R")

respath <- "Results"

#---------------------------------------------------
#  Data reading
#---------------------------------------------------

datatab <- read.table("Data/2_PollutionData.csv", sep = ";", header = T)

# Compute temperature 3-day mean
datatab$T3d <- ma(datatab$Tmean, 3)

# Time variable
datatab$Date <- as.Date(datatab$Date)
datatab$numday <- as.numeric(datatab$Date)

datatab <- na.omit(datatab)
n <- nrow(datatab)

#---------------------------------------------------
#     Applying CGAIM
#---------------------------------------------------

cgaim_res <- cgaim(MCV ~ 
  g(NO2, O3, PM2.5, bs = "mpi", constraints = list(sign.const = 1)) + 
  s(numday) + s(T3d), data = datatab,
  alpha.control = list(norm.type = "sum")
)

p <- ncol(cgaim_res$gfit)
d <- ncol(cgaim_res$indexfit)

#---------------------------------------------------
#    Bootstrapping for uncertainty evaluation
#---------------------------------------------------

cl <- makeCluster(10)
cgaim_CI <- confint(cgaim_res, B = B, l = 7, 
  applyFun = "parLapply", cl = cl)
stopCluster(cl)

save(cgaim_res, cgaim_CI, 
  file = sprintf("%s/3_ Bootstrap_results.RData", respath))
load(sprintf("%s/3_ Bootstrap_results.RData", respath))
  
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

dev.print(png, filename = sprintf("%s/3_Fig5.png", respath), res = 60, 
  width = dev.size()[1], height = dev.size()[2], units = "in")
dev.copy2eps(file = sprintf("%s/3_Fig5.eps", respath))
