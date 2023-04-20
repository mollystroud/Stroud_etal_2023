# SpectralResample.R
# Code to resample AVIRIS-NG reflectance data
# Molly Stroud, 2021-2023

################################################################################
# DESCRIPTION
################################################################################

# The following script contains code for resampling reflectance values from 
# the original AVIRIS-NG hyperspectral data (5 nm) to poorer resolutions


################################################################################
# LOAD IN LIBRARIES
################################################################################
library(stats)
################################################################################



################################################################################
# Description:
# Function to take in a file of reflectance data and output a resampled file of 
# reflectance data
################################################################################
spectralresample <- function(file, bandwidth) {
  file <- read_csv(file)
  bw <- bandwidth
  x <- seq(1, ncol(file))
  y <- as.numeric(file[1, 1:111])
  gaussian <- data.frame()
  for (row in 1:(nrow(file))) { # smooth each row in reflectances file
    processed <- t(as.data.frame(ksmooth(x, y, kernel="normal", 
                                         bandwidth=bandwidth, 
                                         n.points=round(22/(bandwidth/5), digits=0), 
                                         x.points=(seq(1, 111, bandwidth/5)))))
    processed <- processed[-1,]
    gaussian <- rbind(gaussian, processed)
    y <- as.numeric(file[(1 + row), 1:111])
  }
  index <- seq(1, ncol(file), by=bandwidth/5)
  aviris <- file[,index] # get appropriate column names
  spectral <- gaussian
  colnames(spectral) <- colnames(aviris)
  write.csv(spectral, 
            paste0("out/", "spectral.", bandwidht, ",", gsub("([0-9]+).*$", "\\1", file), "m.csv"), 
            row.names=FALSE)
}



