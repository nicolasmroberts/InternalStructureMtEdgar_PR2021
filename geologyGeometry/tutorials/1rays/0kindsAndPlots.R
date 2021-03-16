


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

# In mathematics, a ray is a direction in space. Whereas a line extends 
# infinitely in two opposite directions, a ray starts at a point and extends 
# infinitely in one direction from that point. Various geologic data manifest 
# as rays. This tutorial begins our study of the statistics of ray data. By 
# the way, the standard reference is the book 'Directional Statistics' by 
# Mardia and Jupp (2000).



### IF YOU'VE JUST (RE)STARTED R ###

# If you've just restarted R, then you may need to set your working directory. 
# Your current working directory is named next to the word 'Console' at the 
# top of RStudio's Console pane. It should be the 'geologyGeometry' directory, 
# which contains subdirectories 'data', 'library', etc. If it isn't, then go 
# to the Files pane, navigate to the geologyGeometry directory, click the 
# 'More' menu, and choose 'Set As Working Directory'.

# Then click anywhere on the line of code below, and press RStudio's Run button 
# to execute it. What does it do? It loads the geologyGeometry library into R's 
# memory.
source("library/all.R")



### PALEOMAGNETIC DIRECTIONS ###

# Paleomagnetism is a common source of ray data. Paleomagnetic directions 
# indicate not just the alignment of the magnetic field (a line) but also 
# which way is north (hence a ray). Here is a data set of Titus et al. (in 
# prep) from Iceland. Lower-hemisphere rays are plotted as filled circles. 
# Upper-hemisphere rays are plotted as unfilled circles where their antipodes 
# would plot on the lower hemisphere.
icePmag <- geoDataFromFile("data/icelandOurFlowInSitu29March.csv")
rayEqualAreaPlot(icePmag$direction)

# Notice that there seems to be some polarity reversal in those data. When 
# reversals are present, researchers sometimes ignore the polarities, 
# effectively converting the ray data set into a line data set.



### SOME CRYSTALLOGRAPHIC AXES ###

# Here is a data set consisting of 761 quartz orientations (more specifically, 
# alpha-quartz). They were obtained by electron backscatter diffraction (EBSD) 
# from a single grain of a quartzite sample from the Moine thrust, Scotland 
# (Strine and Wojtal, 2004; Michels et al., 2015).
michelsData <- geoDataFromFile("data/moine_one_grainABCxyz.tsv")

# The a-axes of quartz are rays subject to a 3-fold symmetry. The following 
# plot shows them as 3 * 761 dots.
michelsA1Axes <- lapply(michelsData$rotation, function(r) r[1,])
michelsA2Axes <- lapply(michelsData$rotation, function(r) (oriTrigonalTrapezohedralGroup[[2]] %*% r)[1,])
michelsA3Axes <- lapply(michelsData$rotation, function(r) (oriTrigonalTrapezohedralGroup[[3]] %*% r)[1,])
rayEqualAreaPlotThree(michelsA1Axes, michelsA2Axes, michelsA3Axes, colorA="red", colorB="green", colorC="blue")

# However, the crystallographic symmetry of quartz is such that the c-axis is 
# a line, not a ray. If you really want to see the c-axes, here they are.
michelsCAxes <- lapply(michelsData$rotation, function(r) r[3,])
lineEqualAreaPlot(michelsCAxes)

# Foreshadowing: Divorcing each a-axis triple from its associated c-axis 
# represents a loss of information. In the orientation section of these 
# tutorials, we treat these crystallographic orientations holistically as 
# orientations.



### OTHER EXAMPLES ###

# Any plane has a pole line. But if that plane has a preferred side, then that 
# pole is a ray rather than a line. For example, consider a set of bedding 
# planes. If the younging directions are known, then the poles have a 
# preferred sense and are rays. If the younging is unknown, then the poles are 
# lines.

# Suppose that we have a fault with a known slip direction. The pole to the 
# fault is a line, not a ray. But how do we characterize the direction of slip 
# within the fault plane? It turns out the best descriptor is a 'vorticity' 
# vector, which is a ray. For example, suppose that the fault is vertical. If 
# the slip is dextral, then the vorticity ray points down; if the slip is 
# sinistral, then the vorticity ray points up.



### CONCLUSION ###

# Many geologic data types are treatable as rays.


