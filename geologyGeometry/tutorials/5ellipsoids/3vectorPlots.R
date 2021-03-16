


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

# This tutorial demonstrates how to make plots of 5D or 6D log-ellipsoid 
# vectors. The plots, being 2D or 3D, cannot show all of the five or six 
# degrees of freedom. So these plots, like the Hsu-Nadai and orientation plots 
# before them, are quite imperfect. You should always view your ellipsoidal 
# data in several different styles of plot, to lessen your chances of missing 
# important features of your data set.



### PLOTTING LOG-ELLIPSOID VECTORS ###

# Load a data set of 2,536 spinel grain ellipsoids, obtained by X-ray computed 
# tomography from a hand sample from New Caledonia station KA04-17 (Titus et 
# al., 2011; Chatzaras et al., in prep).
newcal <- geoEllipsoidDataFromAvizoFile("data/newcalSpinelKA04-17.tsv", doNormalize=TRUE)

# Previous tutorials have argued that log-ellipsoid vectors are the best way 
# to compute with ellipsoids. So it might behoove us to plot ellipsoids in 
# this format. There's just one problem: The vectors are 5D or 6D, so we can't 
# visualize them.
newcal$vector

# So we have to resign ourselves to plotting 2D or 3D projections. The 
# following command creates a giant plot, showing 10 different 2D projections 
# of the data set. In particular, the plot in the ith row and jth column shows 
# how the ith and jth components of the vectors relate to each other.
ellPairsPlot(newcal$vector)

# Notice that the third column and row have a special feature: The data set 
# seems to be bimodal (meaning that it has two regions of high density). For 
# example, the following command focuses on the 2nd and 3rd components.
ellVectorPlot(c(2, 3), newcal$vector)

# You can't really see this bimodality in any of our usual plots.
ellHsuNadaiPlot(newcal$logA)
ellEqualAreaPlot(newcal$rotation, newcal$a)
ellEqualVolumePlot(newcal$rotation, newcal$a, simplePoints=TRUE)

# By the way, ellVectorPlot can also show three components at a time. For 
# example, here are the 2nd, 3rd, and 5th components. The bimodality is 
# striking.
ellVectorPlot(c(2, 3, 5), newcal$vector, simplePoints=TRUE)



### PLOTTING PRINCIPAL COMPONENTS ###

# Plotting log-ellipsoid vectors is especially powerful when combined with 
# principal component analysis (PCA). Intuitively, PCA translates and rotates 
# the vectors to give you the best possible view of their dispersion.
pca <- ellPrincipalComponentAnalysis(newcal$vector)
pcs <- lapply(newcal$vector, function(v) as.numeric(t(pca$rotation) %*% (v - pca$center)))
ellPairsPlot(pcs)

# Focus particularly on the first two principal components, which usually give 
# the best possible picture of dispersion. Then contour the data by density. 
# It turns out that the second mode contains few enough points that the 
# contouring algorithm doesn't bother with it.
ellVectorPlot(c(1, 2), pcs)
par(new=TRUE)
contour(kde2d(sapply(pcs, function(v) v[[1]]), sapply(pcs, function(v) v[[2]])))

# So maybe it's better described as a long tail than as a second mode. In any 
# case, subsequent analyses, such as those provided by the R package 
# 'mixtools', can dissect the data set further.



### SYNTHETIC EXAMPLE ###

# In an earlier tutorial, we examined a synthetic example of plane-line pairs, 
# in which there was an outlier that couldn't be seen in the equal-area plot. 
# Now we do a similar example for ellipsoids. Load a synthetic data set of 31 
# volume-normalized ellipsoids.
synthOutlierData <- geoDataFromFile("data/synthEllipsoidOutlier.csv")

# Here are the ellipsoid orientations. Does any datum look like an outlier?
ellEqualVolumePlot(synthOutlierData$rotation, synthOutlierData$a)

# Here are the ellipsoid shapes. Does any datum look like an outlier?
ellHsuNadaiPlot(synthOutlierData$logA)

# Here are all 2D projections of the ellipsoid vectors. Does any datum look 
# like an outlier?
ellPairsPlot(synthOutlierData$vector)

# Here are the same three plots, with half of the data arbitrarily colored 
# blue, half colored green, and the outlier colored red. You see that the 
# outlier's orientation is similar to the blue ones, while its shape is 
# similar to the green ones. If you guessed any outlier at all in the vector 
# plots, then you probably guessed this one.
synthColors <- c(replicate(15, "blue"), replicate(15, "green"), "red")
ellEqualVolumePlot(synthOutlierData$rotation, synthOutlierData$a, colors=synthColors)
ellHsuNadaiPlot(synthOutlierData$logA, colors=synthColors)
ellPairsPlot(synthOutlierData$vector, colors=synthColors)



### CONCLUSION ###

# Don't rely on any one plotting system. View your data in several kinds of 
# plot, to lessen your chances of missing important features of your data set.


