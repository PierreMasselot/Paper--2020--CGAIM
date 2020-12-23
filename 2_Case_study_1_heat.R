#######################################################################
#
#                     Real-world case study
#                            paper CGAIM
#
#######################################################################
library(splines)
library(mgcv)
library(devtools) # Only to install the cgaim package
library(parallel)

# Should be installed from github
# install_github("PierreMasselot/cgaim")
library(cgaim)

source("0_Useful_functions.R")

#---------------------------------------------------
#  Data reading
#---------------------------------------------------

# Read data
dataread <- read.table("Data/1_HeatData.csv", sep = ";", header = T)

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
datamod <- c(datatab[c("Death", "dos", "Year")], laggedX)

#---- Fit the cGAIM
# The gj of Tmin, Tmax and Vp are constrained to be increasing
# The alphas are constrained to be decreasing with lag and all positive
# The norm is the L1 to ensure that each index is a weighted average
system.time(
result <- cgaim(
  Death ~ g(Tmin, fcons = "inc", acons = list(monotone = -1, sign.const = 1)) + 
    g(Tmax, fcons = "inc", acons = list(monotone = -1, sign.const = 1)) +
    g(Vp, fcons = "inc", acons = list(monotone = -1, sign.const = 1)) +
    s(dos, s_opts = list(k = 4)) + s(Year, s_opts = list(k = 3)), 
  data = datamod, alpha.control = list(norm.type = "sum"),
  smooth.control = list(optimizer = "efs", sp = rep(0, 5))
))

#---- Confidence intervals through bootstrap
B <- 1000
n <- nrow(datatab) - 2

# Draw bootstrap samples
datablock <- split(1:n, datatab$Date$year[-(1:2)])
bpool <- 1:length(datablock)
bsamples <- replicate(B, sample(datablock, replace = T), simplify = F)
bsamples <- sapply(bsamples, function (b) rep_len(unlist(b), n))

# Compute CIs. May be done without parallel computing (time-consuming)
set.seed(8)
cl <- makeCluster(10)
resCI <- confint(result, bsamples = bsamples, 
  applyFun = "parLapply", cl = cl, type = "boot.pct")
stopCluster(cl) 

# Save bootstrap results to gain time
save(result, resCI, file = sprintf("%s/2_Bootstrap_results.RData", respath))
load(sprintf("%s/2_Bootstrap_results.RData", respath))

  
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
par(mfcol = c(3, 3), mar = c(5, 0, 0, 2) + .1, oma = c(1, 10, 4, 0), 
  cex.lab = 1.8, xpd = NA, cex.axis = 1.2)
for (j in 1:p){
  # Alphas
  indj <- 1:3 + (j-1) * 3  
  plot(0, 0, ylim = c(0,1), xaxt = "n", col = "white", xlim = c(0.5, 3.5),
    ylab = ifelse(j == 1, expression(alpha), ""),
    yaxt = ifelse(j == 1, "s", "n"), xpd = NA, xlab = "")
  axis(1, at = 1:3, labels = sprintf("Lag %i", 1:3))
  arrows(x0 = 1:3 - .2, y0 = resCI$alpha$boot.bca[indj,1], 
    y1 = resCI$alpha$boot.bca[indj,2], 
    angle = 90, length = 0.05, lwd = 2, code = 3, xpd = T)
  points(1:3 - .2, result$alpha[[j]], pch = 21, bg = cols[1], cex = 3, lwd = 1.5)
  if (j < 3) points(1:3 + .1, classical_alphas, pch = 23, bg = cols[2], cex = 3,
    lwd = 1.5)
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
    ylim = c(0, max(resCI$beta$boot.pct[1:p, 2] * 1.2)), 
    names.arg = "", xpd = NA)
  arrows(x0 = bp, y0 = resCI$beta$boot.bca[j, 1], 
    y1 = resCI$beta$boot.bca[j, 2], 
    angle = 90, length = 0.05, lwd = 2, code = 3, xpd = T)
  
  if (j == 2) legend(par("usr")[1:2], par("usr")[3] - c(0.5, 1.5), 
    c("CGAIM", "", "Classical", ""), pt.bg = rep(cols, each = 2), 
    col = c(1, cols[1], 1, cols[2]), 
    lwd = c(1, 2, 1, 2), lty = c(NA, 1, NA, 2), 
    pch = c(21, NA, 23, NA), bty = "n", ncol = 2, cex = 1.8)
}

dev.print(png, filename = "Results/Figure4.png", res = 200, 
  width = dev.size()[1], height = dev.size()[2], units = "in")
# dev.print(pdf, file = "Results/Figure4.pdf")