#######################################################################
#
#                     Real-world case study
#                            paper CGAIM
#
#######################################################################
library(splines)
library(mgcv)

source("C:/Users/masselpl/Documents/Recherche/2017-2019 - Post-doc/Programmes/R/2 - Indices PPR/gaim/aim.R")
source("C:/Users/masselpl/Documents/Recherche/2017-2019 - Post-doc/Programmes/R/2 - Indices PPR/gaim/Secondary functions.R")

setwd("C:/Users/masselpl/Documents/Recherche/2017-2019 - Post-doc/Programmes/R/2 - Indices PPR/Paper--2020--CGAIM")

#---------------------------------------------------
#  Data reading
#---------------------------------------------------

# Read data
dataread <- read.table("CMMdata.csv", sep = ";", header = T)

# Keep only summer season
summer <- 5:9
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
#   Applying CGAIM
#---------------------------------------------------

result <- aim(y = datatab$Death, 
  x = list(H = datatab[,c("Tmin", "Tmax", "Vp")], dos = datatab$dos, 
    year = datatab$Year),
  smooth.control = list(shape = c("mpi", "tp", "tp")),
  alpha.control = list(norm.type = "L1"),
  trace = FALSE)
datatab$CGAIM <- result$z[,1]




#---------------------------------------------------
#   Computing other Humidity-Temperature variables
#---------------------------------------------------

# Humidex
datatab$Humidex <- datatab$Tmean + (0.5555 * ((datatab$Vp * 10^-2) - 10)) 

# Modelling Humidex
HIres <- gam(Death ~ s(Humidex) + s(dos) + s(Year), data = datatab)

