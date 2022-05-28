library(splines)
library(mgcv)
library(devtools)
library(parallel)
library(forecast)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(doParallel)


# Should be installed from github
# install_github("PierreMasselot/cgaim")
library(cgaim)

source("0_Useful_functions.R")

#---------------------------------------------------
#  Data reading
#---------------------------------------------------

datatab <- read.table("Data/2_PollutionData.csv", sep = ";", header = T)

# Compute temperature 3-day mean
datatab$T3d <- ma(datatab$Tmean, 3)

# Time variable
datatab$Date <- as.Date(datatab$Date)
datatab$dos <- as.numeric(format(datatab$Date, "%j"))
datatab$year <- as.numeric(format(datatab$Date, "%Y"))
datatab$numday <- as.numeric(datatab$Date)

datatab <- na.omit(datatab)
n <- nrow(datatab)

#---------------------------------------------------
#     Applying CGAIM
#---------------------------------------------------

cres <- cgaim(MCV ~ 
    g(NO2, O3, PM2.5, fcons = "inc", acons = list(sign = 1),
      s_opts = list(k = 10), label = "Pol") + 
    s(dos, s_opts = list(k = 7)) + s(year, s_opts = list(k = 3)),
  data = datatab, smooth_control = list(sp = rep(0, 3))
)

ures <- cgaim(MCV ~ 
    g(NO2, O3, PM2.5, label = "Pol") + 
    s(dos, s_opts = list(k = 7)) + s(year, s_opts = list(k = 3)),
  data = datatab, smooth_control = list(sp = rep(0, 3))
)

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

alphadf <- rbind(cdf, udf)
colnames(alphadf)[4:5] <- c("low", "high")
alphadf$pol <- factor(alphadf$pol)

# create data.frame for curves
# cgdf <- data.frame(method = "CGAIM", 
#   index = c(cres$indexfit), 
#   est = c(cres$gfit[,names(cres$alpha)]))
# 
# ugdf <- data.frame(method = "GAIM", 
#   index = c(ures$indexfit), 
#   est = c(ures$gfit[,names(ures$alpha)]))
# 
# gdf <- rbind(cgdf, ugdf)
# 
# cgcidf <- data.frame(method = "CGAIM", cci$g[,1,], cres$indexfit[,1])
# colnames(cgcidf)[2:4] <- c("low", "high", "z")
# 
# ugcidf <- data.frame(method = "GAIM", uci$g[,1,], ures$indexfit[,1])
# colnames(ugcidf)[2:4] <- c("low", "high", "z")
# 
# gcidf <- rbind(cgcidf, ugcidf)

#---- Plot results

# Color palette
pal <- rev(brewer.pal(3, "Blues")[-1])
names(pal) <- c("CGAIM", "GAIM")

# Pollutant labels
pollab <- c(NO2 = expression(NO[2]), O3 = expression(O[3]), 
  "PM2.5" = expression(PM[2.5]))

# X
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


# aplot <- ggplot(alphadf, aes(x = pol, group = method, color = method)) + 
#   theme_classic() + 
#   geom_pointrange(aes(y = est, ymin = low, ymax = high), 
#     position = position_dodge(width = .5), size = .5) +
#   geom_hline(yintercept = 0) + 
#   scale_color_manual(name = "", 
#     values = pal) + 
#   scale_x_discrete(name = "") + 
#   scale_y_continuous(name = expression(alpha)) + 
#   theme(panel.grid.major.y = element_line(colour = "grey", linetype = 2),
#     panel.border = element_rect(colour = grey(.9), fill = NA))

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


# gplot <- ggplot() + theme_classic() +
#   # geom_ribbon(aes(x = z, ymin = low, ymax = high, fill = method),
#   #   data = gcidf, alpha = .3) +
#   geom_line(aes(x = index, y = est, color = method), data = gdf,
#     size = 1) + 
#   geom_line(aes(x = z, y = low, color = method), linetype = 3,
#     data = gcidf, size = .5) +
#   geom_line(aes(x = z, y = high, color = method), linetype = 3,
#     data = gcidf, size = .5) +
#   scale_color_manual(name = "", 
#     values = pal, guide = "none") + 
#   scale_x_continuous(name = "Index") + 
#   scale_y_continuous(name = "g") + 
#   geom_hline(yintercept = 0) + 
#   theme(panel.grid.major.y = element_line(colour = "grey", linetype = 2),
#     panel.border = element_rect(colour = grey(.9), fill = NA))

# Put together and save
# aplot / gplot + plot_layout(guides = "collect")

# Add legend
par(mar = c(rep(0, 4)))
plot.new()
legend("left", legend = names(pal), col = pal, pch = 16, lwd = 2, bty = "n",
  xpd = T)

# Save
dev.print(pdf, "Figures/SupFigure5.pdf", height = 7)















# pct_cgaim <- cgaim_CI$alpha$boot.pct
# 
# x11(width = 10, height = 5)
# par(mfrow = c(1,2), cex.lab = 1.3, cex.axis = 1.2)
#   
# plot(cgaim_res$alpha[[1]], ylim = c(0,1), col = "white",
#   xlab = "Pollutant", ylab = expression(alpha), xaxt = "n",
#   main = bquote("a) Index weights" ~ italic(alpha)))
# axis(1, at = 1:p, 
#   labels = c(expression(NO[2]), expression(O[3]), expression(PM[2.5])))
# arrows(x0 = 1:p, y0 = pct_cgaim[,1], y1 = pct_cgaim[,2],
#   angle = 90, length = 0.05, lwd = 2, code = 3, col = "darkgrey")
# points(1:p, cgaim_res$alpha[[1]], col = "darkgrey", pch = 21, cex = 2, 
#   bg = "cornflowerblue")
#   
# plot(cgaim_res, select = 1, xlab = "Z", ylab = "g()", 
#   col = "cornflowerblue", lwd = 3,
#   main = bquote("b) Function" ~ italic(g())),
#   ci = cgaim_CI, ci.args = list(col = "lightgrey")
# )
# abline(h = 0)
# 
# dev.print(png, filename = "Results/Figure3.png", res = 200, 
#   width = dev.size()[1], height = dev.size()[2], units = "in")
# # dev.print(pdf, file = "Results/Figure5.pdf")
