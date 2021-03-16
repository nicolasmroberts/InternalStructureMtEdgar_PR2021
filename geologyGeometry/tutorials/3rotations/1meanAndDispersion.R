


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

# In elementary statistics, often the first statistics that one computes for a 
# data set are the sample mean and the sample standard deviation. The mean 
# describes the center of the data set. The standard deviation describes how 
# dispersed or 'spread-out' about that mean the data are.

# This tutorial illustrates the analogous statistics for rotations. It turns 
# out that there are several options for quantifying mean and dispersion.



### EXTRINSIC METHODS ###

# Because of a mathematical coincidence (called 'quaternions'), rotations in 
# 3-dimensional space can be identified with lines in 4-dimensional space. 
# This trick lets us apply a bunch of line techniques to rotations. This 
# section gives two examples.

# First, in elementary statistics, the normal distribution plays a pivotal 
# role. In statistics of rotations, the most popular analogue of the normal 
# distribution is the 'matrix Fisher distribution'. What is it? It's the 
# Bingham distribution applied to lines in four dimensions. Here we use it to 
# make a synthetic data set.
synthTrueMean <- rotUniform()
synthTrueK <- diag(c(15, 7, 1))
synthRots <- rotFisher(m=synthTrueMean, k=synthTrueK, n=100)
rotEqualVolumePlot(synthRots)

# Second, we can describe the mean and dispersion of a data set using the same 
# scatter-matrix approach that we used for lines. Now, the four eigenvectors 
# of the scatter matrix produce four rotations. The first one is the mean. The 
# first two indicate the direction of greatest girdling. The first three give 
# a kind of 'best-fit surface' that describes the data. The fourth is the 
# rotation maximally different from the first three.
synthMeanScat <- rotMeanScatter(synthRots)
rotEqualVolumePlotTwo(synthRots, synthMeanScat$rotations, colorB="red")

# The four eigenvalues of the scatter matrix are non-negative and sum to 1. 
# The first one describes how concentrated the data are (about the mean). The 
# second describes how dispersed the data are (along a girdle through the 
# mean). The third describes how dispersed the data are (along that vague 
# surface). The fourth describes how dispersed the data are (away from the 
# surface).
synthMeanScat$values



### DISTRIBUTIONAL METHODS ###

# In elementary statistics, computing the sample mean and standard deviation 
# is tantamount to fitting a normal distribution to the data. In statistics of 
# rotations, the analogous technique is fitting a matrix Fisher distribution 
# to the data.

# If we assume that our rotation data arise from a matrix Fisher distribution, 
# then this function estimates the mean rotation M and concentration matrix K 
# of that distribution. Here's the estimated M compared to the true M that was 
# used to make the synthetic data set in the first place.
synthMLE <- rotFisherMLE(synthRots)
synthMLE$mHat
synthTrueMean

# It's noteworthy that the mean produced by this technique is always identical 
# to the scatter-matrix mean above.

# Here's the estimated K compared to the true K.
synthMLE$kHat
synthTrueK

# The matrix K is a little hard to interpret, but its eigenvalues are a 
# measure of concentation and their inverses are a measure of dispersion.
eigen(synthMLE$kHat, symmetric=TRUE)$values



### INTRINSIC METHODS ###

# The variance of the data {R1, R2, ..., Rn} about a rotation R is 
# proportional to the sum of the squared distances from R to the Ri. The 
# following code defines a function to help you visualize this concept. You 
# can then give that function any trend, plunge, and amount-of-rotation you 
# want, and observe how the variance changes.
varianceVisualizer <- function(rot) {
  rotEqualVolumePlot(synthRots, curves=lapply(synthRots, rotGeodesicPoints, rot))
  rotVariance(synthRots, rot)
}
varianceVisualizer(rotMatrixFromAxisAngle(c(geoCartesianFromTrendPlungeDeg(c(220, 45)), 120 * degree)))

# The Frechet mean is defined to be the R that minimizes the variance.
synthMeanVar <- rotMeanVariance(synthRots)
varianceVisualizer(synthMeanVar$mean)
synthAA <- rotAxisAngleFromMatrix(synthMeanVar$mean)
geoTrendPlungeDegFromCartesian(synthAA[1:3])
synthAA[[4]] / degree

# The projected arithmetic mean (green in plot below) is usually close to, but 
# not identical to, the Frechet mean (blue in plot).
rotEqualVolumePlotThree(synthRots, list(synthMeanScat$rotations[[1]]), list(synthMeanVar$mean), colorB="green", colorC="blue")
rotDistance(synthMeanScat$rotations[[1]], synthMeanVar$mean) / degree

# The variance of the data about their Frechet mean is a scalar measure of 
# dispersion.
synthMeanVar$variance

# Another technique for quantifying dispersion is principal component analysis 
# (PCA). It is appropriate when the data are tightly concentrated about their 
# mean. Intuitively, it fits a girdle to express how the point cloud spreads 
# out from the mean, and a direction perpendicular to the girdle expressing 
# how the data spread out from the girdle, and the direction perpendicular to 
# those two.
synthPCA <- rotLeftPrincipalComponentAnalysis(synthRots, synthMeanVar$mean, numPoints=10)
rotEqualVolumePlot(synthRots, curves=synthPCA$curves, simplePoints=TRUE)

# The PCA offers three 'magnitudes', which quantify how dispersed the data are 
# in those three directions.
synthPCA$magnitudes



### WHEN NOT TO COMPUTE THE MEAN ###

# The mean is a good summary of the 'center' of a data set only when the data 
# set has a meaningful center. For rays and lines, we've seen that the mean is 
# a bad idea for multi-modal data sets. The same is true for rotations.

# But let's return to our Farallon-Pacific plate motion data set (Engebretson 
# et al., 1984; Prentice, 1987) for different kind of example.
source("data/farpacRotations.R")
rotEqualVolumePlot(farpacRs, colors=hues(farpacTs))

# Recall that these rotations are measured cumulatively relative to the 
# present. Therefore the rotations in the distant past (blues and magentas) 
# implicitly incorporate the rotations of the recent past (reds and oranges). 
# So these rotations are in no way 18 independent measurements of a single 
# phenomenon whose mean we'd like to know. Rather, they are highly 
# interdependent on each other and on time. The sample mean is not going to 
# capture any of that. Although it can be computed easily, it should not be 
# reported as a descriptor of the data. (In other tutorials, we use regression 
# to study the time dependence for this data set.)



### CONCLUSION ###

# To describe your data set, you typically want to plot the data, compute the 
# mean, and compute some measure of dispersion. You have various options.


