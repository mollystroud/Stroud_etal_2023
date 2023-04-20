# PLSR_MonteCarlo.R
# Code to build a PLSR model and run a Monte Carlo simulation using resampled
# AVIRIS-NG data paired with in situ TSS samples
# Molly Stroud, 2021-2023

################################################################################
# DESCRIPTION
################################################################################

# The following script contains code that builds a PLSR model using the PRESS
# statistic, then applies this model to resampled AVIRIS-NG data to analyze how
# TSS retrieval accuracy changes with varying resolutions


################################################################################
# LOAD IN LIBRARIES
################################################################################
library(ggplot2)
library(pls)
library(reshape2)
library(qpcR)
library(ggpmisc)
library(wesanderson)
library(rgl)
library(data.table)
library(plotly)
library(patchwork)
################################################################################
source('src/SNRSimulation.R')

################################################################################
# Function:
# Input a file and get lowest PRESS value for PLSR model
################################################################################
aviris <- read_csv('AVIRIS_reflectances.csv')
colnames(aviris) <- sub("\\..*", "", colnames(asd_data))

press_comps <- function(file) {
  data <- read_csv(file) # read in file
  colnames(data) <- sub("\\..*", "", colnames(data))
  t.colnames <- colnames(data[c(2:ncol(data))]) # save column names for indexing
  v.colnames <- colnames(data[1])
  t.matrix <- as.matrix(aviris[t.colnames]) # create test matrix
  v.matrix <- as.matrix(aviris[v.colnames]) # create validation matrix
  myplsr <- plsr(v.matrix ~ t.matrix, data=data, ncomp = 4, scale=FALSE, validation='LOO') # create model
  press <- myplsr$validation[6] # grab PRESS value
  correct_comps <- which.min(unlist(press))
  return(correct_comps) # return optimal # of components
}

################################################################################
# Function:
# Input a file and get error statistics from the PLSR model
################################################################################
error_stats <- data.frame() # empty data frame to put outputs in
plsr_function <- function(filename) {
  file <- read_csv(filename) # open up file
  colnames(file) <- sub("\\..*", "", colnames(file))
  comps <- press_comps(filename) # get appropriate # of components
  t.colnames <- colnames(file[2:ncol(file)]) # save column names for indexing
  v.colnames <- colnames(file[1])
  t.matrix <- as.matrix(aviris[t.colnames]) # create test matrix
  v.matrix <- as.matrix(v[v.colnames]) # create validation matrix
  myplsr <- plsr(v.matrix ~ t.matrix, ncomp = comps, scale=FALSE, validation='LOO') # create model
  # Get matrix of values to test against predicted values
  x_val_matrix <- as.matrix(file[2:ncol(file)])
  y_val_matrix <- as.matrix(file[1])
  # Get predicted values
  preds <- as.matrix(data.frame(predict(myplsr, ncomp=comps, x_val_matrix/pi)))
  # Calculate R2
  r2 <- mean(diag(cor(preds, as.numeric(y_val_matrix)))^2)
  # Calculate RMSE
  rmse <- sqrt(mean((preds - as.numeric(y_val_matrix))^2))
  # Calculate nRMSE
  nrmse <- rmse/mean(as.numeric(y_val_matrix))
  # bind to dataframe
  error_stats <<- rbind(error_stats, c(regmatches(filename, regexpr("[0-9].*[0-9]", filename)), r2, rmse, nrmse))
}

################################################################################
# Generate SNR on original AVIRIS-NG data at varying levels, rerun 1000 times
# and average out stats
################################################################################
aviris_5m <- read_csv("in/reflectances_5m_all.csv") # open up original, unnoised data file
aviris_5m <- aviris_5m[,-1] # remove sample concentrations, only noise reflectances
samples <- read_csv("in/reflectances_5m_all.csv") # keep file with sample concentrations
for (j in 1:1000) { # 1000 runs
  for (i in seq(50, 1000, 10)) { # snr range from 10 to 1000, intervals of 10
    output <- snr_degradation(aviris_5m, i) # add noise
    output <- cbind(samples[,1], output) # bind sample concentrations and noised data
    colnames(output)[1] <- 'Samples'
    write_csv(output, paste0("out/", i, '_snr_reflectances.csv')) # save out noised file
  }
  snr_files <- list.files("out", pattern="reflectances.csv", all.files=TRUE, full.names=TRUE) # get all noised files
  for(file in snr_files) {
    plsr_function(file) # run model on each noised file and append stats to df
  }
}



################################################################################
# Generate SNR on spatially and spectrally resampled AVIRIS-NG data at varying 
# levels, rerun 1000 times and average out stats
################################################################################
resamp_files <- list.files("out/", pattern="nm", all.files=TRUE, full.names=FALSE)

for(k in resamp_files) {
  file <- read_csv(paste0('out/', k))
  for (j in 1:1000) { # 1000 runs
    for (i in seq(50, 1000, 50)) { # snr range from 10 to 1000, intervals of 10
      output <- snr_degradation(file[,-1], i) # add noise
      output <- cbind(samples[,1], output) # bind sample concentrations and noised data
      colnames(output)[1] <- 'Samples'
      write_csv(output, paste0("out/snr/", 'snr.', i, '.', k)) # save out noised file
    }
    snr_files <- list.files("out/snr/", pattern=k, all.files=TRUE, full.names=TRUE) # get all noised files
    for(h in snr_files) {
      plsr_function(h) # run model on each noised file and append stats to df
    }
  }
  colnames(error_stats) <- c('filename', 'r2', 'rmse', 're')
  write_csv(error_stats, paste0('out/snr/', k, '_snr_stats.csv'))
}

################################################################################
# Get calculated statistics and organized data
################################################################################
stats_files <- list.files("out/snr/", pattern="stats.csv", all.files=TRUE, full.names=T)
all_stats <- data.frame()
for(file in stats_files) {
  data <- read_csv(file)
  means <- aggregate(data[,2:4], list(data$filename), mean)# take means of each respective SNR grouping
  colnames(means)[1] <- 'filename'
  means$snr <- sub("\\..*", "", means$filename) # pull out snr
  spectral <- strsplit(means$filename, "[.]") # pull out spectral
  means$spectral <- sapply(spectral, '[[', 3)
  means$spectral <- gsub('.{2}$', '', means$spectral)
  means$spatial <- sub('.*\\.', '', means$filename) # pull out spatial
  all_stats <- rbind(means, all_stats)
}

stats <- aggregate(all_stats[,2:7], list(all_stats$filename), mean)
write_csv(stats, 'out/snr/all_stats.csv')



################################################################################
# Create heatmaps of nRMSE
################################################################################
stats <- read_csv('out/snr/all_stats.csv')
pal <- wes_palette("Zissou1", 100, type = "continuous")
theme_set(theme_classic())

# spatial v. spectral
model <- lm(nrmse ~ I(spatial^2) + log(spectral), data=stats) # model of best fit
# create grid with resolutions of interest 
grd_spatspec <- CJ(A = seq(5, 100, 5), B = seq(5, 100, 5))
colnames(grd_spatspec) <- c('spatial', 'spectral')
# make predictions using model
grd_spatspec$predictions <- predict(model, grd_spatspec)

# plot
nrmse_spatspec <- ggplot(grd_spatspec, aes(x = spatial, y = spectral, fill = predictions)) +
  geom_tile() + 
  scale_fill_gradientn(colours = (pal), name='nRMSE') +
  xlab("Pixel Size (m)") + ylab("Bandwidth (nm)") + 
  scale_x_continuous(breaks = seq(0, 100, by = 10), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 100, by = 10), expand = c(0, 0)) +
  coord_fixed() + 
  theme(legend.position = "top")


# spatial v. snr
stats <- stats[,-1]
stats_spatsnr <- stats[stats$spectral==5,]
model <- lm(nrmse ~ log(snr) + I(spatial^2), data=stats_spatsnr) # model of best fit
# create grid with resolutions of interest
grd_spatsnr <-CJ(A = seq(10, 1000, 10), B = seq(5, 100, 5))
colnames(grd_spatsnr) <- c('snr', 'spatial')
# make predictions
grd_spatsnr$predictions <- predict(model, grd_spatsnr)
nrmse_spatsnr <- ggplot(grd_spatsnr, aes(x = spatial, y = snr, fill = predictions)) +
  geom_tile() + 
  scale_fill_gradientn(colours = (pal), name='nRMSE') +
  xlab("Pixel Size (m)") + ylab("SNR") + 
  scale_x_continuous(breaks = seq(0, 100, by = 10), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 1000, by = 100), expand = c(0, 0)) +
  theme(legend.position = "bottom")

# spectral v. snr
stats_specsnr <- stats[stats$spatial==5,]
model <- lm(nrmse ~ log(snr) + log(spectral), data=stats_specsnr) # model of best fit
# create grid with resolutions of interest
grd_specsnr <-CJ(A = seq(10, 1000, 10), B = seq(5, 100, 5))
colnames(grd_specsnr) <- c('snr', 'spectral')
# make predictions
grd_specsnr$predictions <- predict(model, grd_specsnr)
nrmse_specsnr <- ggplot(grd_specsnr, aes(x = spectral, y = snr, fill = predictions)) +
  geom_tile() + 
  scale_fill_gradientn(colours = (pal), name='nRMSE') +
  xlab("Bandwidth (nm)") + ylab("SNR") + 
  scale_x_continuous(breaks = seq(0, 100, by = 10), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 1000, by = 100), expand = c(0, 0)) +
  theme(legend.position = "bottom")


nrmse_spatspec + nrmse_spatsnr + nrmse_specsnr
