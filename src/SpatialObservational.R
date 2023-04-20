# SpatialObservational.R
# Code to threshold GRWL by observable river area
# Molly Stroud, 2021-2023

################################################################################
# DESCRIPTION
################################################################################

# The following script contains code that uses GRWL (Allen and Pavelsky 2018)
# to threshold river surface area


################################################################################
# LOAD IN LIBRARIES
################################################################################
library(foreign)
library(MASS)
library(ggplot2)
library(ggrepel)


# read in dbf files, combine, and save out
setwd("in/GRWL_vector_V01.01")
dbfs <- list.files(".", pattern='.dbf', full.names = F)
all_files <- data.frame()
for(file in dbfs) {
  data <- as.data.frame(read.dbf(file))
  all_files <- rbind(all_files, data)
}
write_csv(all_files, 'out/GRWL_alldata.csv')


################################################################################
#### Observability in terms of river area (km2)  
################################################################################

# Allen & Pavelsky, 2018 (GRWL):
grwl  = read_csv('out/GRWL_alldata.csv')
# widths
w_grwl = grwl$width[grwl$lakeFlag == 0 & grwl$elev_m > 0 & grwl$width >= 30] # does not exclude braided rivers
w = w_grwl[w_grwl>=90]
# lengths
meanL = mean(c(sqrt(2)*30, 30))
# area
A = w*meanL*1e-6


# set up observational thresholds:
int = 50
minBreak = floor(min(w)/int)*int
maxBreak = ceiling(1000/int)*int
ResBreaks = seq(minBreak, maxBreak, by=int)
AreaBreaks =  (ResBreaks*1e-3)^2

N = length(AreaBreaks)
obsA = rep(NA, N)
for (i in 1:N){
  print(i)
  obsA[i] = sum(A[AreaBreaks[i]<A])
}

# make data plotable
data <- as.data.frame(cbind(ResOfInterest, modObsAofInterest))
data1 <- as.data.frame(cbind(xSeq, modObsA))

# plot
ggplot(data=data, aes(x=ResOfInterest, y=modObsAofInterest)) + 
  geom_line(data=data1, aes(x=xSeq, y=modObsA), size=1, color='#12A1BA') + 
  geom_point(size=2, fill='#12A1BA', color='black', pch=21) +
  geom_text(aes(label=paste0(ResOfInterest, " m, ", round(modObsAofInterest), " km2")), hjust=-.1) +
  labs(x='Resolution (m)', y=expression(Observeable ~ river ~ area ~ (km^2))) +
  theme_bw()

