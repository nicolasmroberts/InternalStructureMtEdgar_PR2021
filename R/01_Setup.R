# This file loads the necessary libraries, functions, and runs some import and setup scripts

#Load local/custom libraries and functions
source("geologyGeometry/library/all.R") # geologyGeometry is written by Joshua R. Davis
source("R/functions_AMS.R") # Written by Nicolas Roberts
source("R/functions_equal area nets.R") # Written by Nicolas Roberts
source("R/functions_equal volume plots.R") # Written by Nicolas Roberts
source("R/functions_stat custom methods.R") # Written by Nicolas Roberts

# Load Rpackages. If you have not already installed these at some point, use the command install.packages("packageName") for each package
library("ggplot2")
library("plotly")
library("sf")
library("gridExtra")
library("ggforce")
library("ggspatial")
library("raster")
library("scales")
library("imager")
library("tidyverse")
library("grid")
library("cowplot")
library("ggsn")
library("ggnewscale")

# Run data import and setup plots
source("R/02_importMaps.R") # Import the shape files and set the map colors
source("R/03_importAMS.R") # Import AMS data
source("R/04_importVSM.R") # Import magnetic hysteresis data
source("R/05_computeStationMeans.R") # Compute the mean AMS ellipsoid for each station
source("R/06_Join AMS and geology.R") # Perform a "Join" of the AMS data and the geologic map shape file. This will append geology information onto each AMS ellipsoid

