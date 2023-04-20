# SNRSimulation.R
# Code to simulate varying SNRs on input reflectance data
# Molly Stroud, 2021-2023

################################################################################
# DESCRIPTION
################################################################################

# The following script contains code for simulating specified SNR levels on
# reflectance values from  AVIRIS-NG hyperspectral data


################################################################################
# Function:
#Simulate given SNR on a single value
################################################################################
simulate_snr <- function(input, SNR) { # input file, desired SNR
  n = rnorm(n=length(input), mean=0, sd=1) # calculate random normal noise
  Es = sum(input^2)
  En = sum(n^2)
  alpha = sqrt(Es/(SNR*En)) # calculate a scalar used to yield a predefined SNR
  output = input+alpha*n # calculate output based on input and calculated noise
}

################################################################################
# Function:
# Add noise to an entire data frame
################################################################################
add_noise <- function(file, snr) {
  df <- data.frame() # create empty data frame
  for (row in 1:nrow(file)) { # for each row, add noise and bind to data frame
    new_row <- simulate_snr(file[row,], snr)
    df <- rbind(df, new_row)
  }
  return(df)
}








