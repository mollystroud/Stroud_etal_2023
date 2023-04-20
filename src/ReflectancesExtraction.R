# ReflectancesExtraction.R
# Code to extract reflectances from spatially resampled AVIRIS-NG data
# Molly Stroud, 2021-2023

################################################################################
# DESCRIPTION
################################################################################

# The following script contains code for extracting reflectance values from 
# multiple AVIRIS-NG hyperspectral data files that have been resampled to
# different spatial resolutions


################################################################################
# LOAD IN LIBRARIES
################################################################################
library(raster)
library(sp)
library(sf)
library(dplyr)
library(rgdal)
library(spatstat)
library(terra)
library(SpatialEpi)
library(stringr)
################################################################################

# Get locations and values of in situ TSS samples from Delta-X mission
latlong <- read.csv("ASD_Spectra_Pts.csv")
latlong <- latlong[,-c(3,4,5,6)]
colnames(latlong) <- c("Latitude", "Longitude", "TSS")

################################################################################
# Description:
# Create function that inputs a tif, latitude, and longitude, and extracts 
# reflectances from the tif and export to new file
################################################################################
extract_values <- function(tif, lat, long, sample) {
  reflectances <- data.frame() # create empty data frame
  point <- cbind(long, lat) # bind lat long
  ras <- (stack(tif)) # create raster stack from file
  reflectances <- rbind(raster::extract(x=ras, y=point, method='simple'), reflectances) # Extract reflectances, bind
  write.csv(reflectances, paste0(tif, sample, ".csv"), row.names = FALSE) # export csv with reflectance data
}

`%notin%` <- Negate(`%in%`) # create 'not in' function

################################################################################
# Hard coded, creates list of files and extracts reflectances from each grouping
# using function
################################################################################
filenames_194554_all <- list.files("MASKED\\17", pattern="194554", all.files=TRUE, full.names=TRUE)
filenames_194554_hdr <- list.files("MASKED\\17", pattern=glob2rx("194554*.hdr"), full.names=TRUE)
filenames_194554 <- filenames_194554_all[(filenames_194554_all %notin% filenames_194554_hdr)]
lapply(filenames_194554, extract_values, lat=29.41803, long=-91.34284, sample="_28_69")

filenames_203526_all <- list.files("MASKED\\17", pattern="203526", all.files=TRUE, full.names=TRUE)
filenames_203526_hdr <- list.files("MASKED\\17", pattern=glob2rx("203526*.hdr"), full.names=TRUE)
filenames_203526 <- filenames_203526_all[(filenames_203526_all %notin% filenames_203526_hdr)]
lapply(filenames_203526, extract_values, lat=29.69770, long=-91.21597, sample="_27_76")
lapply(filenames_203526, extract_values, lat=29.47420, long=-91.48590, sample="_29_55")

filenames_211636_all <- list.files("MASKED\\17", pattern="211636", all.files=TRUE, full.names=TRUE)
filenames_211636_hdr <- list.files("MASKED\\17", pattern=glob2rx("211636*.hdr"), full.names=TRUE)
filenames_211636 <- filenames_211636_all[(filenames_211636_all %notin% filenames_211636_hdr)]
lapply(filenames_211636, extract_values, lat=29.54604, long=-91.42773, sample="_22_94") 
lapply(filenames_211636, extract_values, lat=29.62887, long=-91.39939, sample="_32_45") 
lapply(filenames_211636, extract_values, lat=29.65046, long=-91.39769, sample="_32_93")
lapply(filenames_211636, extract_values, lat=29.64110, long=-91.39450, sample="_35_76")
lapply(filenames_211636, extract_values, lat=29.65821, long=-91.38839, sample="_45_01")


filenames_202455_all <- list.files("MASKED\\17", pattern="202455", all.files=TRUE, full.names=TRUE)
filenames_202455_hdr <- list.files("MASKED\\17", pattern=glob2rx("202455*.hdr"), full.names=TRUE)
filenames_202455 <- filenames_202455_all[(filenames_202455_all %notin% filenames_202455_hdr)]
lapply(filenames_202455, extract_values, lat=29.65137, long=-91.24182, sample="_22_4")
lapply(filenames_202455, extract_values, lat=29.65640, long=-91.24868, sample="_38_4")


filenames_162140_all <- list.files("MASKED\\18", pattern="162140", all.files=TRUE, full.names=TRUE)
filenames_162140_hdr <- list.files("MASKED\\18", pattern=glob2rx("162140*.hdr"), full.names=TRUE)
filenames_162140 <- filenames_162140_all[(filenames_162140_all %notin% filenames_162140_hdr)]
lapply(filenames_162140, extract_values, lat=29.59740, long=-91.37167, sample="_21_75") 
lapply(filenames_162140, extract_values, lat=29.60882, long=-91.38647, sample="_29_96")
lapply(filenames_162140, extract_values, lat=29.52788, long=-91.26845, sample="_42_77") 
lapply(filenames_162140, extract_values, lat=29.59482, long=-91.33984, sample="_50_61") 


filenames_153319_all <- list.files("MASKED\\18", pattern="153319", all.files=TRUE, full.names=TRUE)
filenames_153319_hdr <- list.files("MASKED\\18", pattern=glob2rx("153319*.hdr"), full.names=TRUE)
filenames_153319 <- filenames_153319_all[(filenames_153319_all %notin% filenames_153319_hdr)]
lapply(filenames_153319, extract_values, lat=29.41683, long=-91.34377, sample="_29_47")


filenames_151140_all <- list.files("MASKED\\18", pattern="151140", all.files=TRUE, full.names=TRUE)
filenames_151140_hdr <- list.files("MASKED\\18", pattern=glob2rx("151140*.hdr"), full.names=TRUE)
filenames_151140 <- filenames_151140_all[(filenames_151140_all %notin% filenames_151140_hdr)]
lapply(filenames_151140, extract_values, lat=29.65738, long=-91.24044, sample="_32_65")


filenames_145243_all <- list.files("MASKED\\18", pattern="145243", all.files=TRUE, full.names=TRUE)
filenames_145243_hdr <- list.files("MASKED\\18", pattern=glob2rx("145243*.hdr"), full.names=TRUE)
filenames_145243 <- filenames_145243_all[(filenames_145243_all %notin% filenames_145243_hdr)]
lapply(filenames_145243, extract_values, lat=29.70033, long=-91.21698, sample="_37_44")


filenames_161313_all <- list.files("MASKED\\18", pattern="161313", all.files=TRUE, full.names=TRUE)
filenames_161313_hdr <- list.files("MASKED\\18", pattern=glob2rx("161313*.hdr"), full.names=TRUE)
filenames_161313 <- filenames_161313_all[(filenames_161313_all %notin% filenames_161313_hdr)]
lapply(filenames_161313, extract_values, lat=29.47007, long=-91.27598, sample="_45_82")

filenames_165604_all <- list.files("MASKED\\18", pattern="165604", all.files=TRUE, full.names=TRUE)
filenames_165604_hdr <- list.files("MASKED\\18", pattern=glob2rx("165604*.hdr"), full.names=TRUE)
filenames_165604 <- filenames_165604_all[(filenames_165604_all %notin% filenames_165604_hdr)]
lapply(filenames_165604, extract_values, lat=29.62887, long=-91.39939, sample="_62_99")
################################################################################


################################################################################
# Hard coded, for each spatial resolution, accumulate the reflectance values and
# bind into a single file for later use
################################################################################
filenames_5m_17 <- list.files("MASKED\\17", pattern=glob2rx("*_5m*.csv"), all.files=TRUE, full.names=TRUE)
filenames_5m_18 <- list.files("MASKED\\18", pattern=glob2rx("*_5m*.csv"), all.files=TRUE, full.names=TRUE)
filenames_5m_all <- append(filenames_5m_17, filenames_5m_18)

filenames_10m_17 <- list.files("MASKED\\17", pattern=glob2rx("*_10m*.csv"), all.files=TRUE, full.names=TRUE)
filenames_10m_18 <- list.files("MASKED\\18", pattern=glob2rx("*_10m*.csv"), all.files=TRUE, full.names=TRUE)
filenames_10m_all <- append(filenames_10m_17, filenames_10m_18)

filenames_15m_17 <- list.files("MASKED\\17", pattern=glob2rx("*_15m*.csv"), all.files=TRUE, full.names=TRUE)
filenames_15m_18 <- list.files("MASKED\\18", pattern=glob2rx("*_15m*.csv"), all.files=TRUE, full.names=TRUE)
filenames_15m_all <- append(filenames_15m_17, filenames_15m_18)

filenames_20m_17 <- list.files("MASKED\\17", pattern=glob2rx("*_20m*.csv"), all.files=TRUE, full.names=TRUE)
filenames_20m_18 <- list.files("MASKED\\18", pattern=glob2rx("*_20m*.csv"), all.files=TRUE, full.names=TRUE)
filenames_20m_all <- append(filenames_20m_17, filenames_20m_18)

filenames_25m_17 <- list.files("MASKED\\17", pattern=glob2rx("*_25m*.csv"), all.files=TRUE, full.names=TRUE)
filenames_25m_18 <- list.files("MASKED\\18", pattern=glob2rx("*_25m*.csv"), all.files=TRUE, full.names=TRUE)
filenames_25m_all <- append(filenames_25m_17, filenames_25m_18)

filenames_30m_17 <- list.files("MASKED\\17", pattern=glob2rx("*_30m*.csv"), all.files=TRUE, full.names=TRUE)
filenames_30m_18 <- list.files("MASKED\\18", pattern=glob2rx("*_30m*.csv"), all.files=TRUE, full.names=TRUE)
filenames_30m_all <- append(filenames_30m_17, filenames_30m_18)

filenames_35m_17 <- list.files("MASKED\\17", pattern=glob2rx("*_35m*.csv"), all.files=TRUE, full.names=TRUE)
filenames_35m_18 <- list.files("MASKED\\18", pattern=glob2rx("*_35m*.csv"), all.files=TRUE, full.names=TRUE)
filenames_35m_all <- append(filenames_35m_17, filenames_35m_18)

filenames_40m_17 <- list.files("MASKED\\17", pattern=glob2rx("*_40m*.csv"), all.files=TRUE, full.names=TRUE)
filenames_40m_18 <- list.files("MASKED\\18", pattern=glob2rx("*_40m*.csv"), all.files=TRUE, full.names=TRUE)
filenames_40m_all <- append(filenames_40m_17, filenames_40m_18)

filenames_50m_17 <- list.files("MASKED\\17", pattern=glob2rx("*_50m*.csv"), all.files=TRUE, full.names=TRUE)
filenames_50m_18 <- list.files("MASKED\\18", pattern=glob2rx("*_50m*.csv"), all.files=TRUE, full.names=TRUE)
filenames_50m_all <- append(filenames_50m_17, filenames_50m_18)

filenames_60m_17 <- list.files("MASKED\\17", pattern=glob2rx("*_60m*.csv"), all.files=TRUE, full.names=TRUE)
filenames_60m_18 <- list.files("MASKED\\18", pattern=glob2rx("*_60m*.csv"), all.files=TRUE, full.names=TRUE)
filenames_60m_all <- append(filenames_60m_17, filenames_60m_18)

filenames_70m_17 <- list.files("MASKED\\17", pattern=glob2rx("*_70m*.csv"), all.files=TRUE, full.names=TRUE)
filenames_70m_18 <- list.files("MASKED\\18", pattern=glob2rx("*_70m*.csv"), all.files=TRUE, full.names=TRUE)
filenames_70m_all <- append(filenames_70m_17, filenames_70m_18)

filenames_80m_17 <- list.files("MASKED\\17", pattern=glob2rx("*_80m*.csv"), all.files=TRUE, full.names=TRUE)
filenames_80m_18 <- list.files("MASKED\\18", pattern=glob2rx("*_80m*.csv"), all.files=TRUE, full.names=TRUE)
filenames_80m_all <- append(filenames_80m_17, filenames_80m_18)

filenames_90m_17 <- list.files("MASKED\\17", pattern=glob2rx("*_90m*.csv"), all.files=TRUE, full.names=TRUE)
filenames_90m_18 <- list.files("MASKED\\18", pattern=glob2rx("*_90m*.csv"), all.files=TRUE, full.names=TRUE)
filenames_90m_all <- append(filenames_90m_17, filenames_90m_18)

filenames_100m_17 <- list.files("MASKED\\17", pattern=glob2rx("*_100m*.csv"), all.files=TRUE, full.names=TRUE)
filenames_100m_18 <- list.files("MASKED\\18", pattern=glob2rx("*_100m*.csv"), all.files=TRUE, full.names=TRUE)
filenames_100m_all <- append(filenames_100m_17, filenames_100m_18)
################################################################################


################################################################################
# Description:
# Function to take every file in each spatial resolution group, rename and 
# organize, and save out with corresponding TSS sample concentration
################################################################################
# Get clean column names
columnnames <- read_csv("colnames.csv")
organize_data <- function(list) {
  data <- data.frame()
  for(item in list) {
    contents <- read.csv(item)
    colnames(contents) <- colnames(columnnames)
    contents <- cbind(item, contents)
    data <- rbind(data, contents)
  }
  data$item <- substr(data$item, 22, 26)
  data$item <- gsub('\\.', '', data$item)
  data$item <- sub('_', '.', data$item)
  names(data)[names(data) == 'item'] <- 'Samples'
  write_csv(data, paste0(list, ".csv"))
}
organize_data(reflectances_5m_all)
organize_data(reflectances_10m_all)
# continue for all spatial res of interest

################################################################################