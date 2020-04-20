#######################################################################
#
#                     Real-world case study
#                            paper CGAIM
#
#######################################################################
library(splines)
library(mgcv)
library(devtools)
library(scam)
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

# Read data
dataread <- read.table("CMMdata.csv", sep = ";", header = T)

# Keep only summer season
summer <- 6:8
datatab <- subset(dataread, Month %in% summer)

# Create time variable
datatab$Date <- as.POSIXlt(with(datatab, 
  sprintf("%i-%02.0f-%02.0f", Year, Month, Day)))
# Day of season variable 
datatab$dos <- sequence(tapply(datatab$Date, datatab$Year,length))
# Day of week variable
datatab$dow <- weekdays(as.POSIXlt(datatab$Date))

# Mean temperature variables
datatab$Tmean <- rowMeans(datatab[,c("Tmin", "Tmax")])


#---------------------------------------------------
#        3-day Tmin, Tmax and Vp indices
#---------------------------------------------------

#---- Prepare data
# Scale data for estimation
Xvars <- datatab[,c("Tmin", "Tmax", "Vp")]
p <- length(Xvars)
# Lag data for indices
laggedX <- lapply(Xvars, dlag, 0:2)
for (j in 1:p){
  colnames(laggedX[[j]]) <- sprintf("l%i", 0:2)
}
datamod <- data.frame(datatab, laggedX)

#---- Fit the cGAIM
# The gj of Tmin, Tmax and Vp are constrained to be increasing
# The alphas are constrained to be decreasing with lag and all positive
# The norm is the L1 to ensure that each index is a weighted average
result <- cgaim(Death ~ g(Tmin.l0, Tmin.l1, Tmin.l2, bs = "mpi", 
    constraints = list(monotone = -1, sign.const = 1)) + 
  g(Tmax.l0, Tmax.l1, Tmax.l2, bs = "mpi", 
    constraints = list(monotone = -1, sign.const = 1)) +
  g(Vp.l0, Vp.l1, Vp.l2, bs = "mpi", 
    constraints = list(monotone = -1, sign.const = 1)) +
  dos + Year, data = datamod, algo.control = list(keep.trace = TRUE),
  alpha.control = list(norm.type = "sum", init.type = "random"),
  smooth.control = list(optimizer = "efs")
)

#---- Confidence intervals through bootstrap
B <- 2000
n <- nrow(datamod) - 2

# Draw bootstrap samples
datablock <- split(1:n, datatab$Date$year[-(1:2)])
bpool <- 1:length(datablock)
bsamples <- replicate(B, sample(datablock, replace = T), simplify = F)
bsamples <- sapply(bsamples, function (b) rep_len(unlist(b), n))

# Compute CIs
cl <- makeCluster(2)
resCI <- confint(result, bsamples = bsamples, 
  applyFun = "parLapply", cl = cl)
stopCluster(cl) 

if (F) save(bsamples, bResults, jackResults, 
  file = sprintf("%s/Bootstrap_results_1.RData", respath))
if (F) load(sprintf("%s/Bootstrap_results_1.RData", respath))

  
#---- Fit a GAM on classical indices
classical_alphas <- c(0.4, 0.4, 0.2)
classical_indices <- sapply(laggedX[1:2], "%*%", classical_alphas)

classical_data <- datatab
classical_data[,c("Tmin", "Tmax")] <- classical_indices
classical_fit <- gam(Death ~ s(Tmin) + s(Tmax) + s(dos) + s(Year), 
  data = classical_data)
classical_gz <- scale(predict(classical_fit, type = 'terms')[,1:2])
classical_se <- predict(classical_fit, type = 'terms', se.fit = T)$se.fit[,1:2] 
classical_indices <- na.omit(classical_indices)

#---- Plot results
cols <- c("cornflowerblue", "indianred")

ylims <- range(c(classical_gz - 2 * classical_se, 
  classical_gz + 2 * classical_se, 
  resCI$g$normal[,1:3,]))

x11(title = "App1_results", width = 10, height = 10)
par(mfcol = c(3, 3), mar = c(5, 0, 0, 2) + .1, oma = c(1, 5, 4, 0), 
  cex.lab = 1.8, xpd = NA, cex.axis = 1.2)
for (j in 1:p){
  # Alphas
  indj <- 1:3 + (j-1) * 3
  bp <- barplot(result$alpha[[j]], space = .1,
    ylim = c(0,1), col = cols[1], border = NA, 
    ylab = ifelse(j == 1, expression(alpha), ""),
    names.arg = sprintf("Lag %i", 1:3),
    yaxt = ifelse(j == 1, "s", "n"), xpd = NA)
  arrows(x0 = bp, y0 = resCI$alpha$boot.bca[indj,1], 
    y1 = resCI$alpha$boot.bca[indj,2], 
    angle = 90, length = 0.05, lwd = 2, code = 3, xpd = T)
  if (j < 3) points(bp, classical_alphas, pch = 21, bg = cols[2], cex = 3)
  mtext(colnames(Xvars)[j], xpd = T, line = 2, cex = 1.5)
  
  # Ridge functions
  plot(result, select = j, ci = resCI, type = "l", lwd = 3, 
    xlab = "Z", ylab = ifelse(j == 1, "g", ""), 
    col = cols[1], ylim = ylims, yaxt = ifelse(j == 1, "s", "n"),
    ci.args = list(col = transparency(cols[1], 0.8), border = NA))
  if (j < 3){
    czord <- order(classical_indices[,j])
    lines(classical_indices[czord,j], classical_gz[czord,j], lwd = 3,
      col = cols[2], lty = 2)
    polygon(c(classical_indices[czord,j], rev(classical_indices[czord,j])),
      c(classical_gz[czord,j] + 2 * classical_se[czord,j], 
        rev(classical_gz[czord,j] - 2 * classical_se[czord,j])),
      col = transparency(cols[2], 0.8), border = NA)       
  } 
  
  # Magnitudes beta
  bp <- barplot(result$beta[j + 1], col = "cornflowerblue", 
    cex.names = 1.5, yaxt = ifelse(j == 1, "s", "n"), 
    ylab = ifelse(j == 1, expression(beta), ""),
    ylim = c(0, max(resCI$beta$boot.bca[1:p, 2] * 1.2)), 
    names.arg = "", xpd = NA)
  arrows(x0 = bp, y0 = resCI$beta$boot.bca[j, 1], 
    y1 = resCI$beta$boot.bca[j, 2], 
    angle = 90, length = 0.05, lwd = 2, code = 3, xpd = T)
  
  if (j == 2) legend(par("usr")[1:2], par("usr")[3] - c(0.5, 1.5), 
    c("CGAIM", "Classical"), col = cols, pt.bg = cols,lty = 1:2, lwd = 2, 
    pch = c(NA, 21), bty = "n", ncol = 2, cex = 2)
}

dev.print(png, filename = sprintf("%s/App1_Results.png", respath), res = 100, 
  width = dev.size()[1], height = dev.size()[2], units = "in")