


# Copyright 2016 Joshua R. Davis
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

# We use measures of orientation dispersion to implement a method called 
# crystallographic vorticity axis (CVA) analysis (Michels et al., 2015).



### PROBLEM ###

# If you have just restarted R, then set the working directory to the 
# 'geologyGeometry' directory, which contains subdirectories 'data', 
# 'library', etc. Execute the following line of code to load our custom 
# library and its dependencies.
source("library/all.R")

# Reddy and Buchan (2005) proposed that, as a rock deforms, its mineral grains 
# generally rotate about the vorticity axis of the deformation. Therefore, in 
# a thin section of a deformed rock, one should be able to infer the direction 
# of vorticity from how the crystallographic axes are dispersed. Michels et 
# al. (2015) implemented this idea in an objective, quantitative way, using 
# orientation statistics.

# To see how, let's return to this data set of 761 alpha-quartz orientations 
# (6-fold symmetry) from a single quartzite grain from the Moine thrust zone, 
# obtained by EBSD (Strine and Wojtal, 2004; Michels et al., 2015). In an 
# equal-area plot, the 010 (green) and 001 (blue) axes seem to be dispersed 
# along small circles with a common pole that is nearly vertical. The 100 
# (red) axes are not as clearly dispersed along a small circle, but that could 
# be because they are close to the pole.
michelsData <- geoDataFromFile("data/moine_one_grainABCxyz.tsv")
michels100s <- lapply(michelsData$rotation, function(r) r[1,])
michels010s <- lapply(michelsData$rotation, function(r) r[2,])
michels001s <- lapply(michelsData$rotation, function(r) r[3,])
lineEqualAreaPlotThree(michels100s, michels010s, michels001s, colorA="red", colorB="green", colorC="blue")

# So the question is: How do you fit circles to all three sets of axes 
# simultaneously? To someone familiar with rotation statistics, the obvious 
# answer is to bind the axes together into orientations and fit a curve to 
# those orientations. So, for each alpha-quartz orientation, we put its three 
# axes into the rows of a matrix, which is then subject to trigonal 
# trapezohedral symmetry.
oriEqualAnglePlot(michelsData$rotation, group=oriTrigonalTrapezohedralGroup, simplePoints=TRUE)



### SOLUTION ###

# The orientations are tightly concentrated, so it is reasonable to work with 
# one symmetric copy and ignore the others. Compute the mean orientation and 
# choose representative rotations near it. Then the principal components 
# capture the main directions of dispersion. In particular, the first 
# principal component produces a geodesic that best-fits the cloud of 
# orientations. Zoom in on this plot to see the geodesics.
michelsMean <- oriFrechetMean(michelsData$rotation, group=oriTrigonalTrapezohedralGroup)
michelsRots <- oriNearestRepresentatives(michelsData$rotation, michelsMean, group=oriTrigonalTrapezohedralGroup)
michelsPCA <- rotLeftPrincipalComponentAnalysis(michelsRots, michelsMean, numPoints=5)
oriEqualAnglePlot(michelsRots, curves=michelsPCA$curves, group=oriTrigonalTrapezohedralGroup, simplePoints=TRUE)

# Now we dissect the first principal geodesic, along its rows, into three 
# curves of directions. They are small circles about the first principal 
# direction.
michelsFirst <- lapply(michelsPCA$curves[[1]], function(r) r[1,])
michelsSecond <- lapply(michelsPCA$curves[[1]], function(r) r[2,])
michelsThird <- lapply(michelsPCA$curves[[1]], function(r) r[3,])
n <- nrow(michelsData)
geoTrendPlungeDegFromCartesian(michelsPCA$directions[,1])
lineEqualAreaPlot(points=c(michels100s, michels010s, michels001s, list(michelsPCA$directions[,1])),
                  colors=c(replicate(n, "red"), replicate(n, "green"), replicate(n, "blue"), "black"),
                  curves=list(michelsFirst, michelsSecond, michelsThird))

# This grain is one of 74,724 grains in a single sample. Each grain produces a 
# dispersion axis in this manner. See the figure in the notes. So the sample 
# has 74,724 dispersion axes. They end up being widely dispersed, but the mean 
# is at trend-plunge (228, 84) degrees. So Michels et al. (2015) took that 
# direction as the overall vorticity direction of the deformation that 
# deformed this sample. And this direction was similar to estimates of 
# vorticity direction from earlier studies.

# Michels et al. (2015) analyzed two other samples, from the Gem Lake shear 
# zone and the western Idaho shear zone. They found that their technique 
# reproduced earlier estimates of vorticity axis well in all three cases.



### CONCLUSION ###

# Crystallographic vorticity axis analysis, when reduced to mathematical/
# statistical machinery, is just principal component analysis.


