library(splines)
library(mgcv)
library(devtools)
library(scam)
library(forecast)
library(parallel)
library(vioplot)

setwd("C:/Users/masselpl/Documents/Recherche/2017-2019 - Post-doc/Programmes/R/2 - Indices PPR/Paper--2020--CGAIM")

load_all("C:/Users/masselpl/Documents/Recherche/2017-2019 - Post-doc/Programmes/R/2 - Indices PPR/gaim")
source("Useful_functions.R")

respath <- "C:/Users/masselpl/Documents/Recherche/2017-2019 - Post-doc/Resultats/Part 2 - GAIM/Article V3"

#---------------------------------------------------
#  Data reading
#---------------------------------------------------

datatab <- read.table("Data2 - Pollution.csv", sep = ";", header = T)

# Compute temperature 3-day mean
datatab$T3d <- ma(datatab$Tmean, 3)

datatab <- na.omit(datatab)
n <- nrow(datatab)

#---------------------------------------------------
#     Applying CGAIM
#---------------------------------------------------
y <- datatab$MCV
x <- data.matrix(datatab[,c("NO2", "O3", "PM2.5", "Date", "T3d")])

gaim_res <- gaim_gn(y = y, x = x, 
  w = rep(1/n, n), index = rep(1,3),
  alpha.control = list(norm.type = "sum"))
gaim_eps <- y - gaim_res$fitted
  
cgaim_res <- gaim_gn(y = y, x = x, 
  w = rep(1/n, n), index = rep(1,3),
  alpha.control = list(norm.type = "sum", Cmat = diag(3)),
  smooth.control = list(shape = c("mpi", "tp", "tp")))
cgaim_eps <- y - cgaim_res$fitted

p <- ncol(gaim_res$gz)
d <- length(gaim_res$index)

#---------------------------------------------------
#    Bootstrapping for uncertainty evaluation
#---------------------------------------------------

#datatab$Year <- substr(datatab$Date, 1, 4)

# Draw bootstrap samples
B <- 2000
#bsamples <- replicate(B, sample(unique(datatab$Year), replace = T))
#bsamples <- matrix(sample(seq_len(n), n * B, replace = T), nrow = n, ncol = B)
l <- 7
bsamples <- sample(seq_len(n - l), (n * B) / l, replace = T)
bsamples <- matrix(sapply(bsamples, function(b) b:(b + l)), nrow = n, ncol = B)
btable <- apply(bsamples, 2, function(b) table(c(b, seq_len(n))) - 1)

# Initialize cluster for parallel computation
cl <- makeCluster(2)
clusterExport(cl, ls())
clusterEvalQ(cl, {
  library(splines)
  library(mgcv)
  library(devtools)
  library(scam)
  library(forecast)
  load_all("C:/Users/masselpl/Documents/Recherche/2017-2019 - Post-doc/Programmes/R/2 - Indices PPR/gaim")
  source("Useful_functions.R")
})

# Apply CGAIM on each bootstrap sample
boot_gaim <- parApply(cl, bsamples, 2, function(b){
  beps <- gaim_eps[b]
  bres <- gaim_gn(y = gaim_res$fitted + beps, x = x, 
    w = rep(1/n, n), index = rep(1,3),
    alpha.control = list(norm.type = "sum"))
  return(bres[c("alpha", "beta", "z", "gz")])
})

boot_cgaim <- parApply(cl, bsamples, 2, function(b){
  beps <- gaim_eps[b]
  bres <- gaim_gn(y = cgaim_res$fitted + beps, x = x, 
    w = rep(1/bn, bn), index = rep(1,3),
    alpha.control = list(norm.type = "sum", Cmat = diag(3)),
    smooth.control = list(shape = c("mpi", "tp", "tp")))
  return(bres[c("alpha", "beta", "z", "gz")])
})

#--- jackknife calculation of acceleration for BCA
mr <- 5
m <- 20
r <- n %% m

jack_gaim <- parSapply(cl, seq_len(mr), function(k){
  Imat <- sapply(seq_len_m, sample.int, n = n, size = n - r)
  Iout <- setdiff(seq_len(n), Imat)
  alphamat <- matrix(0, m, d)
  for (j in seq_len_m) {
    Ij <- setdiff(seq_len_m, j)
    ij <- c(c(Imat[Ij, ], Iout))
    alphamat[j,] <- gaim_gn(y = y[ij], x = x[ij,], 
      w = rep(1/length(ij), length(ij)), index = rep(1,3),
      alpha.control = list(norm.type = "sum"))$alpha
  }
  aa <- apply(alphamat, 2, function(u){
    t. <- (mean(u) - u) * (m - 1)
    (1/6) * sum(t.^3)/(sum(t.^2))^1.5
  })
  aa
})

jack_cgaim <- parSapply(cl, seq_len(mr), function(k){
  Imat <- sapply(seq_len_m, sample.int, n = n, size = n - r)
  Iout <- setdiff(seq_len(n), Imat)
  alphamat <- matrix(0, m, d)
  for (j in seq_len_m) {
    Ij <- setdiff(seq_len_m, j)
    ij <- c(c(Imat[Ij, ], Iout))
    alphamat[j,] <- gaim_gn(y = y[ij], x = x[ij,], 
      w = rep(1/length(ij), length(ij)), index = rep(1,3),
      alpha.control = list(norm.type = "sum", Cmat = diag(3)),
    smooth.control = list(shape = c("mpi", "tp", "tp")))$alpha
  }
  aa <- apply(alphamat, 2, function(u){
    t. <- (mean(u) - u) * (m - 1)
    (1/6) * sum(t.^3)/(sum(t.^2))^1.5
  })
  aa
})

stopCluster(cl)

if (F) save(bsamples, boot_gaim, boot_cgaim, jack_gaim, jack_cgaim, 
  file = sprintf("%s/Bootstrap_results_2.RData", respath))
if (F) load(sprintf("%s/Bootstrap_results_2.RData", respath))

  
#---------------------------------------------------
#     Visualizing results
#---------------------------------------------------

#alpha_cis_u <- t(confint_gn(gaim_res, parm = 1:3))
#alpha_cis_c <- t(confint_gn(cgaim_res, parm = 1:3))

#------ Coefficients Alpha -----
# Unconstrained
balphas_u <- sapply(boot_gaim,"[[", "alpha")                     
alpha_means_u <- apply(balphas_u, 1, mean)
alpha_cis_u <- apply(balphas_u, 1, quantile, c(0.025, 0.975))

# Constrained
balphas_c <- sapply(boot_cgaim,"[[", "alpha")
alpha_means_c <- apply(balphas_c, 1, mean)
alpha_cis_c <- apply(balphas_c, 1, quantile, c(0.025, 0.975))

#------ BCA --------

#Blist <- list(Y = t(btable), tt = balphas_u[1,], t0 = gaim_res$alpha[1])
#bca_u <- bcajack2(B = Blist)

accel_gaim <- apply(jack_gaim, 1, mean)
accel_cgaim <- apply(jack_cgaim, 1, mean)

z0_gaim <- mapply(function(tt, t0){qnorm(mean(tt < t0))},
  as.data.frame(t(balphas_u)), gaim_res$alpha)
z0_cgaim <- mapply(function(tt, t0){qnorm(mean(tt < t0))},
  as.data.frame(t(balphas_c)), cgaim_res$alpha)

bca_gaim <- sapply(c(0.025, 0.975), function(level){
  za <- qnorm(level)
  correc <- z0_gaim + (z0_gaim + za) / (1 - accel_gaim * (z0_gaim + za))
  indB <- trunc(pnorm(correc) * B)
  indB <- pmin(pmax(indB, 1), B)
  lim <- mapply(function(x, i) sort(x)[i], as.data.frame(t(balphas_u)), indB)
})

bca_cgaim <- sapply(c(0.025, 0.975), function(level){
  za <- qnorm(level)
  correc <- z0_cgaim + (z0_cgaim + za) / (1 - accel_cgaim * (z0_cgaim + za))
  indB <- trunc(pnorm(correc) * B)
  indB <- pmin(pmax(indB, 1), B)
  lim <- mapply(function(x, i) sort(x)[i], as.data.frame(t(balphas_c)), indB)
})



x11(width = 10, height = 5)
par(mfrow = c(1,2), cex.lab = 1.3, cex.axis = 1.2)
  
plot(cgaim_res$alpha, ylim = c(0,1), col = "white",
  xlab = "Pollutant", ylab = expression(alpha), xaxt = "n",
  main = bquote("a) Index weights" ~ italic(alpha)))
axis(1, at = 1:p, 
  labels = c(expression(NO[2]), expression(O[3]), expression(PM[2.5])))
arrows(x0 = 1:p, y0 = bca_cgaim[,1], y1 = bca_cgaim[,2],
  angle = 90, length = 0.05, lwd = 2, code = 3, col = "darkgrey")
points(1:p, cgaim_res$alpha, col = "darkgrey", pch = 21, cex = 2, 
  bg = "cornflowerblue")
  
plot(cgaim_res$am.fit, select = 1, rug = F, xlab = "Z", ylab = "g()", 
  col = "cornflowerblue", shade = T, lwd = 3,
  main = bquote("b) Function" ~ italic(g())))

dev.print(png, filename = sprintf("%s/App2_Results.png", respath), res = 60, 
  width = dev.size()[1], height = dev.size()[2], units = "in")


x11()
par(mfrow = c(1,2))

bp <- barplot(gaim_res$alpha, ylim = c(0,1), cex.names = 1.3, cex.axis = 1.2,
  names.arg = c(expression(NO[2]), expression(O[3]), expression(PM[2.5])))
arrows(x0 = bp, y0 = alpha_cis_u[1,], y1 = alpha_cis_u[2,],
  angle = 90, length = 0.05, lwd = 2, code = 3)
points(bp, alpha_means_u, cex = 1.5, pch = 16)

bp <- barplot(cgaim_res$alpha, ylim = c(0,1), cex.names = 1.3, cex.axis = 1.2,
  names.arg = c(expression(NO[2]), expression(O[3]), expression(PM[2.5])))
arrows(x0 = bp, y0 = alpha_cis_c[1,], y1 = alpha_cis_c[2,],
  angle = 90, length = 0.05, lwd = 2, code = 3)
points(bp, alpha_means_c, cex = 1.5, pch = 16)

x11()
par(mfrow = c(1,2))

bp <- barplot(gaim_res$alpha, ylim = c(0,1), cex.names = 1.3, cex.axis = 1.2,
  names.arg = c(expression(NO[2]), expression(O[3]), expression(PM[2.5])))
arrows(x0 = bp, y0 = bca_gaim[,1], y1 = bca_gaim[,2],
  angle = 90, length = 0.05, lwd = 2, code = 3)

bp <- barplot(cgaim_res$alpha, ylim = c(0,1), cex.names = 1.3, cex.axis = 1.2,
  names.arg = c(expression(NO[2]), expression(O[3]), expression(PM[2.5])))
arrows(x0 = bp, y0 = bca_cgaim[,1], y1 = bca_cgaim[,2],
  angle = 90, length = 0.05, lwd = 2, code = 3)
points(bp, alpha_means_c, cex = 1.5, pch = 16)

x11()
par(mfrow = c(1,2))

vioplot(as.data.frame(t(balphas_u)), col = "lightgrey", drawRect = F, 
  cex.names = 1.3, cex.axis = 1.2,
  names = c(expression(NO[2]), expression(O[3]), expression(PM[2.5])))
segments(x0 = 1:p - .2, x1 = 1:p + .2, y0 = gaim_res$alpha, lwd = 3)
points(gaim_res$alpha, cex = 2, pch = 16)
points(alpha_means_u, cex = 2, pch = 0)

vioplot(as.data.frame(t(balphas_c)), col = "lightgrey", drawRect = F, 
  cex.names = 1.3, cex.axis = 1.2,
  names = c(expression(NO[2]), expression(O[3]), expression(PM[2.5])))
segments(x0 = 1:p - .2, x1 = 1:p + .2, y0 = cgaim_res$alpha, lwd = 3)
points(cgaim_res$alpha, cex = 2, pch = 16)
points(alpha_means_c, cex = 2, pch = 0)


#------ Ridge functions gz -----
# Unconstrained
bgs_u <- sapply(boot_gaim, function(bres){
  gznorm <- spline(x = fscaling(bres$z), y = bres$gz[,1],
    xout = fscaling(1:1000))$y
})
g_cis_u <- t(apply(bgs_u, 1, quantile, c(0.025, 0.975)))

# Constrained
bgs_c <- sapply(boot_cgaim, function(bres){
  gznorm <- spline(x = fscaling(bres$z), y = bres$gz[,1],
    xout = fscaling(1:1000))$y
})
g_cis_c <- t(apply(bgs_c, 1, quantile, c(0.025, 0.975)))

x11()
par(mfrow = c(1,2))

plot(sort(fscaling(gaim_res$z)), gaim_res$gz[order(gaim_res$z),1], 
  type = "l", lwd = 2,
  xlab = "Index value", ylab = "g(.)", cex.lab = 1.3, cex.axis = 1.2)
matlines(fscaling(1:1000), g_cis_u, lty = 2, col = 1)

plot(sort(fscaling(cgaim_res$z)), cgaim_res$gz[order(cgaim_res$z),1], 
  type = "l", lwd = 2,
  xlab = "Index value", ylab = "g(.)", cex.lab = 1.3, cex.axis = 1.2)
matlines(fscaling(1:1000), g_cis_c, lty = 2, col = 1)