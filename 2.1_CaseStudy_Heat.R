################################################################################
#
#  R code for the real-world application of 
#
#   Masselot et al., 2022
#   Constrained groupwise additive index models
#   Biostatistics
#
#   Section 5 Application
#
#   Author: Pierre Masselot
#
################################################################################

#-------------------------------------------
# Packages
#-------------------------------------------

#----- Used packages
library(RColorBrewer) # Colorpalette Blues
library(tsModel) # For Lag function
library(doParallel) # For parallel

#----- cgaim package
library(cgaim)

#---------------------------------------------------
#  Data reading
#---------------------------------------------------

# Read data
dataread <- read.table("Data/2.1_HeatData.csv", sep = ";", header = T)

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

# Mean temperature variable
datatab$Tmean <- rowMeans(datatab[,c("Tmin", "Tmax")])

#---------------------------------------------------
#        3-day Tmin, Tmax and Vp indices
#---------------------------------------------------

#---- Prepare data

# Select variables
Xvars <- datatab[,c("Tmin", "Tmax", "Vp")]
p <- length(Xvars)

# Lag data for indices
laggedX <- lapply(Xvars, Lag, 0:2)
for (j in 1:p){
  colnames(laggedX[[j]]) <- sprintf("l%i", 0:2)
}

# Final dataset
datamod <- c(datatab[c("Death", "dos", "Year")], laggedX)

# Basic model formula
formbasis <- "Death ~ s(dos, s_opts = list(k = 4)) + s(Year, s_opts = list(k = 3))"

# Submodels to compare
modlist <- expand.grid(0:1, 0:1, 0:1)[-1,]

#---- Fit the cGAIM
# The gj of Tmin, Tmax and Vp are constrained to be increasing
# The alphas are constrained to be decreasing with lag and all positive
# The norm is the L1 to ensure that each index is a weighted average

# Index terms
termlist <- c("g(Tmin, fcons = 'inc', acons = list(monotone = -1, sign = 1))",
  "g(Tmax, fcons = 'inc', acons = list(monotone = -1, sign = 1))", 
  "g(Vp, fcons = 'inc', acons = list(monotone = -1, sign = 1))")

# Fit all submodels
creslist <- apply(modlist, 1, function(i){
  formula <- paste(c(formbasis, termlist[as.logical(i)]), collapse = " + ")
  cgaim(as.formula(formula), data = datamod, 
    control = list(sm_pars = list(sp = rep(0, 2 + sum(i)))))
})

# Keep model with lowest GCV
gcvlist <- sapply(creslist, "[[", "gcv")
cres <- creslist[[6]]

#---- Fit the GAIM

# Index terms
utermlist <- c("g(Tmin)", "g(Tmax)", "g(Vp)")

# Fit all submodels
ureslist <- apply(modlist, 1, function(i){
  formula <- paste(c(formbasis, utermlist[as.logical(i)]), collapse = " + ")
  cgaim(as.formula(formula), data = datamod, 
    control = list(sm_pars = list(sp = rep(0, 2 + sum(i)))))
})

# Keep model with lowest GCV
ugcvlist <- sapply(ureslist, "[[", "gcv")
ures <- ureslist[[6]]

#---- Confidence intervals through bootstrap

# Number of bootstrap samples
B <- 1000

# Confidence intervals for CGAIM
cci <- confint(cres, type = "boot", B = B, nc = detectCores() - 2)

# Confidence intervals for Unconstrained model
uci <- confint(ures, type = "boot", B = B, nc = detectCores() - 2)

#----- Prepare results for plotting

# Create data.frame for alphas
cdf <- data.frame(method = "CGAIM", 
  var = rep(names(cres$alpha), lengths(cres$alpha)), 
  lag = c(sapply(lengths(cres$alpha), seq_len)) - 1, 
  est = unlist(cres$alpha), cci$alpha)
udf <- data.frame(method = "GAIM", 
  var = rep(names(ures$alpha), lengths(ures$alpha)), 
  lag = c(sapply(lengths(ures$alpha), seq_len)) - 1, 
  est = unlist(ures$alpha), uci$alpha)

# Put together
alphadf <- rbind(cdf, udf)
colnames(alphadf)[5:6] <- c("low", "high")

# Change x to put models side by side (manual dodge)
alphadf$x <- alphadf$lag + ifelse(alphadf$method == "CGAIM", -.1, .1)

# List of selected indices
selinds <- unique(alphadf$var)
ninds <- length(selinds)

#----- Figure 5: Results

# Palette
pal <- brewer.pal(3, "Blues")[-1]
names(pal) <- c("GAIM", "CGAIM")

# Layout matrix to control the position of nice labels
laymat <- matrix(0, nrow = 4, ncol = ninds + 1)
laymat[1, -1] <- 1:ninds
laymat[3, -1] <- ninds + 1:ninds
laymat[c(1,3), 1] <- 2*ninds + 1:2
laymat[2, -1] <- 2 * ninds + 3
laymat[4, -1] <- 2 * ninds + 4
laymat <- cbind(laymat, 9)

# Initialize plot
layout(laymat, width = c(.1, rep(1, ninds), .2), heights = rep(c(1, .2), 2))
defpar <- par()
par(mar = c(3, 3, 4, 2))

#----- Plot alphas
for (ind in selinds){
  # Select variable
  inddf <- subset(alphadf, var == ind)
  
  # Plot points
  plot(est ~ x, inddf, col = pal[inddf$method], pch = 16, cex = 1.5,
    ylim = range(alphadf[,c("low", "high")]), xaxt = "n", xlab = "", ylab = "")
  
  # Add grid
  abline(h = axTicks(2), lty = 2, col = "grey")
  abline(h = 0)
  
  # Add segments
  arrows(x0 = inddf$x, y0 = inddf$low, y1 = inddf$high, lwd = 2,
    col = pal[inddf$method], angle = 90, length = .05, code = 3)

  # Add axis
  axis(1, at = unique(alphadf$lag))
  
  # Add title
  title(main = ind)
}

#----- Plot g functions
par(mar = c(3, 3, 1, 2))
for (ind in selinds){
  
  # Compute limits
  ylims <- range(cbind(uci$g[,ind,], cci$g[,ind,]))
  
  # Plot for GAIM
  plot(ures, ci = uci, select = ind, col = pal["GAIM"], lwd = 2,
    ci.plot = "lines", ci.args = list(lty = 2, col = pal["GAIM"]),
    xlab = "Standardized Index", ylim = ylims, ylab = "",
    xcenter = min(ures$indexfit[,ind]),
    xscale = diff(range(ures$indexfit[,ind])))
  
  # Add a grid
  grid()
  
  # Add for CGAIM
  plot(cres, ci = cci, select = ind, col = pal["CGAIM"], lwd = 2,
    ci.plot = "lines", ci.args = list(lty = 2, col = pal["CGAIM"]), add = T, 
    xcenter = min(cres$indexfit[,ind]),
    xscale = diff(range(cres$indexfit[,ind])))
}

#----- Add labels

# remove margins
par(mar = rep(0, 4))

# Add y axis labels
plot.new()
text(par("usr")[2], mean(par("usr")[3:4]), expression(alpha), 
  srt = 90, pos = 2, cex = 1.5, xpd = NA)
plot.new()
text(par("usr")[2], mean(par("usr")[3:4]), "g(.)", 
  srt = 90, pos = 2, cex = 1.5, xpd = NA)

# Add x axis labels
plot.new()
text(mean(par("usr")[1:2]), par("usr")[4], "Lag", 
  pos = 1, cex = 1.5, xpd = NA)
plot.new()
text(mean(par("usr")[1:2]), par("usr")[4], "Standardized Index", 
  pos = 1, cex = 1.5, xpd = NA)

# Add legend
par(mar = c(rep(0, 4)))
plot.new()
legend("left", legend = names(pal), col = pal, pch = 16, lwd = 2, bty = "n",
  xpd = T)
