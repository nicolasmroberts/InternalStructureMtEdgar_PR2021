


# Copyright 2017 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not 
# use this file except in compliance with the License. You may obtain a copy 
# of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software 
# distributed under the License is distributed on an "AS IS" BASIS, WITHOUT 
# WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the 
# License for the specific language governing permissions and limitations 
# under the License.



### INTRODUCTION ###

# A line is tantamount to a pair of rays pointing in opposite directions. 
# Various geologic data types manifest as lines. This tutorial begins our 
# study of line data. By the way, the standard reference is the book 
# 'Directional Statistics' by Mardia and Jupp (2000).



### IF YOU'VE JUST (RE)STARTED R ###

# If you've just restarted R, then you may need to set your working directory. 
# Your current working directory is named next to the word 'Console' at the 
# top of RStudio's Console pane. It should be the 'geologyGeometry' directory, 
# which contains subdirectories 'data', 'library', etc. If it isn't, then go 
# to the Files pane, navigate to the geologyGeometry directory, click the 
# 'More' menu, and choose 'Set As Working Directory'.

# Then click anywhere on the line of code below, and press RStudio's Run button 
# to execute it. It loads the geologyGeometry library into R's memory.
source("library/all.R")



### POLES TO PLANES (DIKES, FOLIATIONS, BEDDING, ETC.) ###

# Load some dike data from Cyprus (Titus, pers. comm.). Make an equal-area, 
# lower-hemisphere plot of the dike poles. The plot appears in the Plots pane. 
# To get a better look, zoom it, resize the window, etc.
site230 <- geoDataFromFile("data/cyprusDikesSite230.tsv")
lineEqualAreaPlot(site230$pole)

# The data file itself contains easting, northing, strike, and dip. The 
# geoDataFromFile function automatically converts the strike and dip into a 
# pole line represented by a unit vector in Cartesian coordinates. Most of 
# the time you don't need to know these details, though.
site230



### LINEATIONS ###

# Load a data set of foliation-lineation pairs from the western Idaho shear 
# zone (Giorgis and Tikoff, 2004).
wiszData <- geoDataFromFile("data/wiszFollins.tsv")

# The lineations are lines (drawn as squares) and the poles to foliations are 
# lines (drawn as circles).
lineEqualAreaPlotTwo(wiszData$direction, wiszData$pole, shapeA="s", shapeB="c")

# Foreshadowing: Because each lineation is paired with a foliation, to treat 
# them separately, as in the preceding plot, is to suffer a loss of 
# information. In the orientation section of these tutorials, we learn how to 
# treat them together.



### SOME CRYSTALLOGRAPHIC AXES ###

# Here is a data set consisting of 761 alpha-quartz orientations. They were 
# obtained by electron backscatter diffraction (EBSD) from a single grain of a 
# quartzite sample. from the Moine thrust, Scotland (Strine and Wojtal, 2004; 
# Michels et al., 2015).
michelsData <- geoDataFromFile("data/moine_one_grainABCxyz.tsv")

# The crystallographic symmetry of quartz is such that the c-axis is a line.
michelsCAxes <- lapply(michelsData$rotation, function(r) r[3,])
lineEqualAreaPlot(michelsCAxes)

# However, the a-axes of quartz are not lines. They are rays subject to a 3-
# fold symmetry. If you really want to see them, here they are.
michelsA1Axes <- lapply(michelsData$rotation, function(r) r[1,])
michelsA2Axes <- lapply(michelsData$rotation, function(r) (oriTrigonalTrapezohedralGroup[[2]] %*% r)[1,])
michelsA3Axes <- lapply(michelsData$rotation, function(r) (oriTrigonalTrapezohedralGroup[[3]] %*% r)[1,])
rayEqualAreaPlotThree(michelsA1Axes, michelsA2Axes, michelsA3Axes, colorA="red", colorB="green", colorC="blue")

# Foreshadowing: Like foliation-lineation pairs, these crystallographic 
# orientations should really be treated holistically as orientations.



### MANY OTHER EXAMPLES ###

# Principal stress directions are lines.

# The directions of an ellipsoid's axes are lines.

# Paleomagnetic directions are rays, not lines. However, when magnetic 
# reversals are present in a data set, it is common to discard the polarity 
# information, effectively converting the rays into lines.



### CONCLUSION ###

# Many kinds of geologic data manifest as lines.


