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

load_all("C:/Users/masselpl/Documents/Recherche/2017-2019 - Post-doc/Programmes/R/2 - Indices PPR/gaim")
source("Useful_functions.R")

respath <- "C:/Users/masselpl/Documents/Recherche/2017-2019 - Post-doc/Resultats/Part 2 - GAIM/Article V3"

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
scdat <- as.data.frame(scale(Xvars))
# Lag data for indices
laggedX <- lapply(scdat, dlag, 0:2, na.action = na.omit)
# Add year and day-of-season to predictors
X <- c(laggedX, list(dos = datatab$dos[-(1:2)], year = datatab$Year[-(1:2)]))
Xall <- Reduce(cbind, X)
# Response
Y <- datatab$Death[-(1:2)]
n <- length(Y)

#---- Fit the cGAIM
# The gj of Tmin, Tmax and Vp are constrained to be increasing
# The alphas are constrained to be decreasing with lag and all positive
# The norm is the L1 to ensure that each index is a weighted average
result <- gaim_gn(y = Y, x = Xall, index = rep(1:3, each = 3), 
  w = rep(1/n, n), 
  smooth.control = list(shape = c("mpi", "mpi", "mpi", "tp", "tp")),
  alpha.control = list(norm.type = "sum", 
    Cmat = const.matrix(index = rep(1:3, each = 3), monotone = rep(-1, 3), 
    sign.const = rep(1, 3))),
  tol = 5e-3, convergence_criterion = "LS", keep.trace = F
)
# ~ 4.7 min
eps <- Y - result$fitted

#---- Confidence intervals through bootstrap
B <- 2000

# Draw bootstrap samples
#bsamples <- matrix(sample(seq_len(n), n * B, replace = T), nrow = n, ncol = B)
datablock <- split(1:n, datatab$Date$year)
bpool <- 1:length(datablock)

# Draw bootstrap samples
bsamples <- replicate(B, sample(datablock, replace = T), simplify = F)
bsamples <- sapply(bsamples, function (b) rep_len(unlist(b), n))
btable <- apply(bsamples, 2, function(b) table(c(b, seq_len(n))) - 1)

# Initialize cluster for parallel computation
cl <- makeCluster(2)
# Transfer objects in cluster
clusterExport(cl, ls())
clusterEvalQ(cl, {
  library(splines)
  library(mgcv)
  library(devtools)
  library(scam)
  library(parallel)
  load_all("C:/Users/masselpl/Documents/Recherche/2017-2019 - Post-doc/Programmes/R/2 - Indices PPR/gaim")
})

# Reapply CGAIM on each bootstrap sample
bResults <- parApply(cl, bsamples, 2, function(b){
  epsb <- eps[b]
  resb <- gaim_gn(y = result$fitted + epsb, x = Xall, 
    index = rep(1:3, each = 3), 
    w = rep(1/n, n), 
    smooth.control = list(shape = c("mpi", "mpi", "mpi", "tp", "tp")),
    alpha.control = list(norm.type = "sum", 
      Cmat = const.matrix(index = rep(1:3, each = 3), monotone = rep(-1, 3), 
      sign.const = rep(1, 3))),
    tol = 5e-3, convergence_criterion = "LS"
  )
  return(resb[c("alpha", "z", "gz", "beta")])
})

# Jackknife for acceleration
mr <- 5
m <- 20
seq_len_m <- seq_len(m)
r <- n %% m
clusterExport(cl, ls())

jackResults <- parSapply(cl, seq_len(mr), function(k){
  Imat <- sapply(seq_len_m, sample.int, n = n, size = n - r)
  Iout <- setdiff(seq_len(n), Imat)
  alphamat <- matrix(0, m, p * 3)
  betamat <- matrix(0, m, p)
  for (j in seq_len_m) {
    Ij <- setdiff(seq_len_m, j)
    ij <- c(c(Imat[Ij, ], Iout))
    mres <- gaim_gn(y = Y[ij], x = Xall[ij,], 
      w = rep(1/length(ij), length(ij)), index = rep(1:3, each = 3),
      smooth.control = list(shape = c("mpi", "mpi", "mpi", "tp", "tp")),
      alpha.control = list(norm.type = "sum", 
          Cmat = const.matrix(index = rep(1:3, each = 3), monotone = rep(-1, 3), 
                                               sign.const = rep(1, 3))),
      tol = 5e-3, convergence_criterion = "LS")
    alphamat[j,] <- mres$alpha
    betamat[j,] <- mres$beta[1:p + 1]
  }
  aa <- apply(alphamat, 2, function(u){
    t. <- (mean(u) - u) * (m - 1)
    (1/6) * sum(t.^3)/(sum(t.^2))^1.5
  })
  ab <- apply(betamat, 2, function(u){
    t. <- (mean(u) - u) * (m - 1)
    (1/6) * sum(t.^3)/(sum(t.^2))^1.5
  })
  c(aa, ab)
})

stopCluster(cl) 

if (F) save(bsamples, bResults, jackResults, 
  file = sprintf("%s/Bootstrap_results_1.RData", respath))
if (F) load(sprintf("%s/Bootstrap_results_1.RData", respath))

bAlphas <- sapply(bResults, "[[", "alpha")
# Remove cases with QP fail
Cout <- apply(bAlphas, 2, function(x) any(x < 0 | x > 1))
bAlphas <- bAlphas[,!Cout]
B <- B - sum(Cout)
#bmeans <- apply(bAlphas, 1, mean)
#cis <- apply(bAlphas, 1, quantile, c(0.025, 0.975))

accel <- apply(jackResults[1:(p*3),], 1, mean)
z0 <- mapply(function(tt, t0){qnorm(mean(tt <= t0))},
  as.data.frame(t(bAlphas)), result$alpha)

bca_ci <- sapply(c(0.025, 0.975), function(level){
  za <- qnorm(level)
  correc <- z0 + (z0 + za) / (1 - accel * (z0 + za))
  indB <- trunc(pnorm(correc) * B)
  indB <- pmin(pmax(indB, 1), B)
  lim <- mapply(function(x, i) sort(x)[i], as.data.frame(t(bAlphas)), indB)
})

#-----  Betas -----

bBetas <- sapply(bResults, "[[", "beta")
baccel <- apply(jackResults[1:3 + p*3,], 1, mean)
bz0 <- mapply(function(tt, t0){qnorm(mean(tt <= t0))},
  as.data.frame(t(bBetas)), result$beta)[1:p + 1]

beta_ci <- sapply(c(0.025, 0.975), function(level){
  za <- qnorm(level)
  correc <- bz0 + (bz0 + za) / (1 - baccel * (bz0 + za))
  indB <- trunc(pnorm(correc) * B)
  indB <- pmin(pmax(indB, 1), B)
  lim <- mapply(function(x, i) sort(x)[i], as.data.frame(t(bBetas))[,1:p + 1], indB)
})
  
#---- Fit a GAM on classical indices
classical_alphas <- c(0.4, 0.4, 0.2)
uns_lagX <- Map(function(x, m, s) x * s + m, laggedX,
  apply(Xvars, 2, mean), apply(Xvars, 2, sd))
classical_indices <- sapply(uns_lagX[1:2], "%*%", classical_alphas)

classical_fit <- gam(Y ~ s(Tmin) + s(Tmax) + s(dos) + s(year), 
  data = data.frame(Y, classical_indices, dos = datatab$dos[-(1:2)], 
    year = datatab$Year[-(1:2)]))
classical_gz <- scale(predict(classical_fit, type = 'terms')[,1:2])
classical_se <- predict(classical_fit, type = 'terms', se.fit = T)$se.fit[,1:2] 

#---- Plot results
cols <- c("cornflowerblue", "indianred")

indices <- mapply("%*%", uns_lagX, split(result$alpha, rep(1:p, each = 3)))
ise <- predict(result$am.fit, type = "terms", se.fit = TRUE)$se.fit[,1:p]

ylims <- range(c(classical_gz - 2 * classical_se, 
  classical_gz + 2 * classical_se, 
  result$gz[,1:p] + 2 * ise, result$gz[,1:p] - 2 * ise))

x11(title = "App1_results", width = 10, height = 10)
par(mfcol = c(3, 3), mar = c(5, 0, 0, 2) + .1, oma = c(1, 5, 4, 0), 
  cex.lab = 1.8, xpd = NA, cex.axis = 1.2)
for (j in 1:p){
  # Alphas
  indj <- 1:3 + (j-1) * 3
  bp <- barplot(result$alpha[indj], space = .1,
    ylim = c(0,1), col = cols[1], border = NA, 
    ylab = ifelse(j == 1, expression(alpha), ""),
    names.arg = sprintf("Lag %i", 1:3),
    yaxt = ifelse(j == 1, "s", "n"), xpd = NA)
  if(j < 3) points(bp, classical_alphas, pch = 21, bg = cols[2], cex = 3)
  arrows(x0 = bp, y0 = bca_ci[indj,1], y1 = bca_ci[indj,2], 
    angle = 90, length = 0.05, lwd = 2, code = 3, xpd = T)
  mtext(colnames(Xvars)[j], xpd = T, line = 2, cex = 1.5)
  
  # Ridge functions
  zord <- order(indices[,j])
  plot(indices[zord,j], result$gz[zord,j], type = "l", lwd = 3, 
    xlab = "Z", ylab = ifelse(j == 1, "g", ""), 
    col = cols[1], ylim = ylims, yaxt = ifelse(j == 1, "s", "n"))
  polygon(c(indices[zord,j], rev(indices[zord,j])),
    c(result$gz[zord,j] + 2 * ise[zord, j], 
      rev(result$gz[zord,j] - 2 * ise[zord, j])),
    col = transparency(cols[1], 0.8), border = NA
  )
  if (j %in% 1:2){
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
    ylim = c(0, max(beta_ci[1:p, 2] * 1.2)), 
    names.arg = "", xpd = NA)
  arrows(x0 = bp, y0 = beta_ci[j, 1], y1 = beta_ci[j, 2], 
    angle = 90, length = 0.05, lwd = 2, code = 3, xpd = T)
  
  if (j == 2) legend(par("usr")[1:2], par("usr")[3] - c(0.5, 1.5), 
    c("CGAIM", "Classical"), col = cols, pt.bg = cols,lty = 1:2, lwd = 2, 
    pch = c(NA, 21), bty = "n", ncol = 2, cex = 2)
}

dev.print(png, filename = sprintf("%s/App1_Results.png", respath), res = 100, 
  width = dev.size()[1], height = dev.size()[2], units = "in")