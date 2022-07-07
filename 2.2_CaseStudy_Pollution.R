################################################################################
#
#  R code for the additional real-world application of 
#
#   Masselot et al., 2022
#   Constrained groupwise additive index models
#   Biostatistics
#
#   Supplementary materials section 2
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

datatab <- read.table("Data/2.2_PollutionData.csv", sep = ";", header = T)

# Time variables
datatab$Date <- as.Date(datatab$Date)
datatab$numday <- as.numeric(datatab$Date)

# Day-of-season variable
datatab$dos <- as.numeric(format(datatab$Date, "%j"))

# Year variable
datatab$year <- as.numeric(format(datatab$Date, "%Y"))

# Remove NAs and keep sample size
datatab <- na.omit(datatab)
n <- nrow(datatab)

#---------------------------------------------------
#     Applying CGAIM
#---------------------------------------------------

# Fit CGAIM
cres <- cgaim(MCV ~ 
    g(NO2, O3, PM2.5, fcons = "inc", acons = list(sign = 1),
      s_opts = list(k = 10), label = "Pol") + 
    s(dos, s_opts = list(k = 7)) + s(year, s_opts = list(k = 3)),
  data = datatab, control = list(sm_pars = list(sp = rep(0, 3)))
)

# Fit unconstrained model
ures <- cgaim(MCV ~ 
    g(NO2, O3, PM2.5, label = "Pol") + 
    s(dos, s_opts = list(k = 7)) + s(year, s_opts = list(k = 3)),
  data = datatab, control = list(sm_pars = list(sp = rep(0, 3)))
)

# Keep number of indices and terms
p <- ncol(cres$gfit)
d <- ncol(cres$indexfit)

#---------------------------------------------------
#    Bootstrapping for uncertainty evaluation
#---------------------------------------------------

# Number of resampling
B <- 1000

# Constrained model
cci <- confint(cres, type = "boot", B = B, nc = max(1, detectCores() - 2))

# Unconstrained model
uci <- confint(ures, type = "boot", B = B, nc = max(1, detectCores() - 2))

#---------------------------------------------------
#     Visualizing results
#---------------------------------------------------

# Create result data.frame
cdf <- data.frame(method = "CGAIM", 
  pol = names(cres$alpha[[1]]), 
  est = unlist(cres$alpha), cci$alpha)
udf <- data.frame(method = "GAIM", 
  pol = names(cres$alpha[[1]]), 
  est = unlist(ures$alpha), 
  uci$alpha)

# Put together and rename
alphadf <- rbind(cdf, udf)
colnames(alphadf)[4:5] <- c("low", "high")
alphadf$pol <- factor(alphadf$pol)

#---- Supplementary Figure 5

# Color palette
pal <- rev(brewer.pal(3, "Blues")[-1])
names(pal) <- c("CGAIM", "GAIM")

# Pollutant labels
pollab <- c(NO2 = expression(NO[2]), O3 = expression(O[3]), 
  "PM2.5" = expression(PM[2.5]))

# Dodge two models
alphadf$x <- as.numeric(alphadf$pol) + 
  ifelse(alphadf$method == "CGAIM", -.1, .1)

# Prepare layout
laymat <- cbind(1:2, 3)
layout(laymat, widths = c(1, .2))
par(mar = c(5, 4, 2, 2))

# Alphas
plot(est ~ x, alphadf, col = pal[alphadf$method], pch = 16, cex = 1.5,
  ylim = range(alphadf[,c("low", "high")]), xaxt = "n", xlab = "Pollutant", 
  ylab = expression(alpha))
abline(h = axTicks(2), lty = 2, col = "grey")
abline(h = 0)
arrows(x0 = alphadf$x, y0 = alphadf$low, y1 = alphadf$high, 
  lwd = 2, col = pal[alphadf$method], angle = 90, length = .05, code = 3)
axis(1, at = seq_along(unique(alphadf$pol)), 
  labels = pollab[levels(alphadf$pol)])

# Ridge function g
ylims <- range(cbind(uci$g[,1,], cci$g[,1,]))
plot(ures, ci = uci, select = 1, col = pal["GAIM"], lwd = 2,
  ci.plot = "lines", ci.args = list(lty = 2, col = pal["GAIM"]),
  xlab = "Standardized index", ylim = ylims, ylab = "",
  xcenter = min(ures$indexfit[,1]),
  xscale = diff(range(ures$indexfit[,1])) / 100)
grid()
plot(cres, ci = cci, select = 1, col = pal["CGAIM"], lwd = 2,
  ci.plot = "lines", ci.args = list(lty = 2, col = pal["CGAIM"]), add = T, 
  xcenter = min(cres$indexfit[,1]),
  xscale = diff(range(cres$indexfit[,1])) / 100)

# Add legend
par(mar = c(rep(0, 4)))
plot.new()
legend("left", legend = names(pal), col = pal, pch = 16, lwd = 2, bty = "n",
  xpd = T)
