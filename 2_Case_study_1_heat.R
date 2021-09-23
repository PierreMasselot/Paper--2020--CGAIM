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
library(ggplot2)
library(RColorBrewer)
library(patchwork)

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

# For models
formbasis <- "Death ~ s(dos, s_opts = list(k = 4)) + s(Year, s_opts = list(k = 3))"
modlist <- expand.grid(0:1, 0:1, 0:1)[-1,]

#---- Fit the cGAIM
# The gj of Tmin, Tmax and Vp are constrained to be increasing
# The alphas are constrained to be decreasing with lag and all positive
# The norm is the L1 to ensure that each index is a weighted average

# Terms
termlist <- c("g(Tmin, fcons = 'inc', acons = list(monotone = -1, sign.const = 1))",
  "g(Tmax, fcons = 'inc', acons = list(monotone = -1, sign.const = 1))", 
  "g(Vp, fcons = 'inc', acons = list(monotone = -1, sign.const = 1))")

creslist <- apply(modlist, 1, function(i){
  formula <- paste(c(formbasis, termlist[as.logical(i)]), collapse = " + ")
  cgaim(as.formula(formula), data = datamod, 
    smooth_control = list(sp = rep(0, 2 + sum(i))))
})
gcvlist <- sapply(creslist, "[[", "gcv")
cres <- creslist[[which.min(gcvlist)]]

# cres <- cgaim(
#   Death ~ g(Tmin, fcons = "inc", acons = list(monotone = -1, sign.const = 1)) + 
#     g(Tmax, fcons = "inc", acons = list(monotone = -1, sign.const = 1)) +
#     g(Vp, fcons = "inc", acons = list(monotone = -1, sign.const = 1)) +
#     s(dos, s_opts = list(k = 4)) + s(Year, s_opts = list(k = 3)), 
#   data = datamod, smooth_control = list(sp = rep(0, 5))
# )

#---- Fit the GAIM
utermlist <- c("g(Tmin)", "g(Tmax)", "g(Vp)")

ureslist <- apply(modlist, 1, function(i){
  formula <- paste(c(formbasis, utermlist[as.logical(i)]), collapse = " + ")
  cgaim(as.formula(formula), data = datamod, 
    smooth_control = list(sp = rep(0, 2 + sum(i))))
})
ugcvlist <- sapply(ureslist, "[[", "gcv")
ures <- ureslist[[which.min(ugcvlist)]]

#---- Confidence intervals through bootstrap
B <- 1000
n <- nrow(datatab) - 2

# Draw bootstrap samples
datablock <- split(1:n, datatab$Date$year[-(1:2)])
bpool <- 1:length(datablock)
bsamples <- replicate(B, sample(datablock, replace = T), simplify = F)
bsamples <- sapply(bsamples, function (b) rep_len(unlist(b), n))

# Compute CIs. May be done without parallel computing (time-consuming)
cl <- makeCluster(max(1, detectCores() - 2))
cci <- confint(cres, bsamples = bsamples, 
  applyFun = "parLapply", cl = cl, type = "boot.pct")
stopCluster(cl) 

# The same for unconstrained model
cl <- makeCluster(max(1, detectCores() - 2))
uci <- confint(ures, bsamples = bsamples, 
  applyFun = "parLapply", cl = cl, type = "boot.pct")
stopCluster(cl) 

# Save bootstrap results to gain time
save(cres, cci, ures, uci, file = "Results/2_Bootstrap_results.RData")
load("Results/2_Bootstrap_results.RData")

  
#---- Fit a GAM on classical indices
# classical_alphas <- c(0.4, 0.4, 0.2)
# classical_indices <- sapply(laggedX[1:2], "%*%", classical_alphas)
# 
# classical_data <- datatab
# classical_data[,c("Tmin", "Tmax")] <- classical_indices
# classical_fit <- gam(Death ~ s(Tmin) + s(Tmax) + s(dos) + s(Year), 
#   data = classical_data)
# classical_gz <- scale(predict(classical_fit, type = 'terms')[,1:2])
# classical_se <- predict(classical_fit, type = 'terms', se.fit = T)$se.fit[,1:2] 
# classical_indices <- na.omit(classical_indices)

# Create result data.frame
cdf <- data.frame(method = "CGAIM", 
  var = rep(names(cres$alpha), lengths(cres$alpha)), 
  lag = c(sapply(lengths(cres$alpha), seq_len)) - 1, 
  est = unlist(cres$alpha), cci$alpha$boot.pct)
udf <- data.frame(method = "GAIM", 
  var = rep(names(ures$alpha), lengths(ures$alpha)), 
  lag = c(sapply(lengths(ures$alpha), seq_len)) - 1, 
  est = unlist(ures$alpha), 
  uci$alpha$boot.pct)

alphadf <- rbind(cdf, udf)
colnames(alphadf)[5:6] <- c("low", "high")

# create data.frame for curves
cgdf <- data.frame(method = "CGAIM", 
  var = rep(names(cres$alpha), each = n), 
  index = c(cres$indexfit), 
  est = c(cres$gfit[,names(cres$alpha)]))

ugdf <- data.frame(method = "GAIM", 
  var = rep(names(ures$alpha), each = n), 
  index = c(ures$indexfit), 
  est = c(ures$gfit[,names(ures$alpha)]))

gdf <- rbind(cgdf, ugdf)

cgcidf <- data.frame(method = "CGAIM", 
  var = rep(names(cres$alpha), each = n),
  apply(cci$g$boot.pct$g[,seq_along(cres$alpha),,drop = F], 3, c),
  c(cci$g$boot.pct$z[,seq_along(cres$alpha)]))
colnames(cgcidf)[3:5] <- c("low", "high", "z")

ugcidf <- data.frame(method = "GAIM", 
  var = rep(names(ures$alpha), each = n), 
  apply(uci$g$boot.pct$g[,seq_along(ures$alpha),,drop = F], 3, c),
  c(uci$g$boot.pct$z[,seq_along(ures$alpha)]))
colnames(ugcidf)[3:5] <- c("low", "high", "z")

gcidf <- rbind(cgcidf, ugcidf)

#---- Plot results

# Alphas
aplot <- ggplot(alphadf, aes(x = lag, group = method, color = method)) + theme_classic() + 
  geom_pointrange(aes(y = est, ymin = low, ymax = high), 
    position = position_dodge(width = .5), size = .5) +
  geom_hline(yintercept = 0) + 
  scale_color_manual(name = "", 
    values = rev(brewer.pal(3, "Blues")[-1])) + 
  scale_x_continuous(name = "Lag", breaks = unique(alphadf$lag)) + 
  scale_y_continuous(name = expression(alpha)) + 
  facet_wrap(~ var) +
  theme(panel.grid.major.y = element_line(colour = "grey", linetype = 2),
    panel.border = element_rect(colour = grey(.9), fill = NA))

# Ridge function g
gplot <- ggplot() + theme_classic() +
  # geom_ribbon(aes(x = z, ymin = low, ymax = high, fill = method),
  #   data = gcidf, alpha = .3) +
  geom_line(aes(x = index, y = est, color = method), data = gdf,
    size = 1) + 
  geom_line(aes(x = z, y = low, color = method), linetype = 3,
    data = gcidf, size = .5) +
  geom_line(aes(x = z, y = high, color = method), linetype = 3,
    data = gcidf, size = .5) +
  scale_color_manual(name = "", 
    values = rev(brewer.pal(3, "Blues")[-1]), guide = "none") + 
  scale_x_continuous(name = "Index") + 
  scale_y_continuous(name = "g") + 
  geom_hline(yintercept = 0) + 
  facet_wrap(~ var, scales = "free") + 
  theme(panel.grid.major.y = element_line(colour = "grey", linetype = 2),
    panel.border = element_rect(colour = grey(.9), fill = NA))

# Put together and save
aplot / gplot + plot_layout(guides = "collect")

ggsave("Figures/Figure4.pdf", height = 7)

# #---- Plot results
# cols <- c("cornflowerblue", "indianred")
# 
# # ylims <- range(c(classical_gz - 2 * classical_se, 
# #   classical_gz + 2 * classical_se, 
# #   resCI$g$boot.pct$g[,1:3,]))
# ylims <- c(-10, 30)
# 
# x11(title = "App1_results", width = 10, height = 10)
# par(mfcol = c(3, 3), mar = c(5, 0, 0, 2) + .1, oma = c(1, 10, 4, 0), 
#   cex.lab = 1.8, cex.axis = 1.2)
# for (j in 1:p){
#   # Alphas
#   indj <- 1:3 + (j-1) * 3  
#   plot(0, 0, ylim = c(0,1), xaxt = "n", col = "white", xlim = c(0.5, 3.5),
#     ylab = ifelse(j == 1, expression(alpha), ""),
#     yaxt = ifelse(j == 1, "s", "n"), xpd = NA, xlab = "")
#   axis(1, at = 1:3, labels = sprintf("Lag %i", 1:3))
#   arrows(x0 = 1:3 - .2, y0 = resCI$alpha$boot.pct[indj,1], 
#     y1 = resCI$alpha$boot.pct[indj,2], 
#     angle = 90, length = 0.05, lwd = 2, code = 3, xpd = T)
#   points(1:3 - .2, result$alpha[[j]], pch = 21, bg = cols[1], cex = 3, lwd = 1.5)
#   if (j < 3) points(1:3 + .1, classical_alphas, pch = 23, bg = cols[2], cex = 3,
#     lwd = 1.5)
#   mtext(colnames(Xvars)[j], xpd = T, line = 2, cex = 1.5)
#   
#   # Ridge functions
#   plot(result, select = j, ci = resCI, type = "l", lwd = 3, 
#     xlab = "Z", ylab = ifelse(j == 1, "g", ""), 
#     col = cols[1], ylim = ylims, yaxt = ifelse(j == 1, "s", "n"),
#     ci.args = list(col = transparency(cols[1], 0.8), border = NA))
#   if (j < 3){
#     czord <- order(classical_indices[,j])
#     lines(classical_indices[czord,j], classical_gz[czord,j], lwd = 3,
#       col = cols[2], lty = 2)
#     polygon(c(classical_indices[czord,j], rev(classical_indices[czord,j])),
#       c(classical_gz[czord,j] + 2 * classical_se[czord,j], 
#         rev(classical_gz[czord,j] - 2 * classical_se[czord,j])),
#       col = transparency(cols[2], 0.8), border = NA)       
#   }
#   abline(h = 0)
#   
#   # Magnitudes beta
#   bp <- barplot(result$beta[j + 1], col = "cornflowerblue", 
#     cex.names = 1.5, yaxt = ifelse(j == 1, "s", "n"), 
#     ylab = ifelse(j == 1, expression(beta), ""),
#     ylim = c(0, max(resCI$beta$boot.pct[1:p, 2] * 1.2)), 
#     names.arg = "", xpd = NA)
#   arrows(x0 = bp, y0 = resCI$beta$boot.pct[j, 1], 
#     y1 = resCI$beta$boot.pct[j, 2], 
#     angle = 90, length = 0.05, lwd = 2, code = 3, xpd = T)
#   
#   if (j == 2) legend(par("usr")[1:2], par("usr")[3] - c(0.5, 1.5), 
#     c("CGAIM", "", "Classical", ""), pt.bg = rep(cols, each = 2), 
#     col = c(1, cols[1], 1, cols[2]), 
#     lwd = c(1, 2, 1, 2), lty = c(NA, 1, NA, 2), 
#     pch = c(21, NA, 23, NA), bty = "n", ncol = 2, cex = 1.8)
# }
# 
# dev.print(png, filename = "Results/Figure2.png", res = 200, 
#   width = dev.size()[1], height = dev.size()[2], units = "in")
# # dev.print(pdf, file = "Results/Figure4.pdf")