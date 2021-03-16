


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

# This tutorial explains how to compute the mean and measure the dispersion of 
# a data set of ellipsoids. Good news: It's easy.



### MEAN ###

# Load the Cyprus AMS data set again (Titus et al., in prep).
cyprusAMS <- geoEllipsoidDataFromIRMFile("data/cyprus_AMS_groupF.tsv", doNormalize=FALSE)

# Compute the mean by converting to log-ellipsoid vectors, taking the mean, 
# and converting back.
cyprusMean <- ellEllipsoidFromVector(arithmeticMean(cyprusAMS$vector))
cyprusMean

# It's that easy. But we've packaged the calculation into a single convenient 
# function anyway.
cyprusMean <- ellMean(cyprusAMS$vector)
cyprusMean

# Here is the mean in a few plots.
ellHsuNadaiPlot(c(cyprusAMS$logA, list(cyprusMean$logA)), es=0.03,
                colors=c(replicate(length(cyprusAMS$logA), "black"), "red"))
ellEqualAreaPlot(c(cyprusAMS$rotation, list(cyprusMean$rotation)),
                 c(cyprusAMS$logA, list(cyprusMean$logA)),
                 colors=c(replicate(length(cyprusAMS$logA), "black"), "red"))
ellEqualVolumePlot(c(cyprusAMS$rotation, list(cyprusMean$rotation)),
                   c(cyprusAMS$logA, list(cyprusMean$logA)),
                   colors=c(replicate(length(cyprusAMS$logA), "white"), "red"))

# A question for you to ponder: In the Hsu-Nadai plot, the mean plots closer 
# to the vertex than do most of the ellipsoids that it averages. What does 
# this mean, geometrically? Why should it be true?



### NICE PROPERTIES IN SPECIAL CASES ###

# If all of your ellipsoids have the same volume (perhaps because they've been 
# volume-normalized), then the mean has that same volume.
cyprusAMS <- geoEllipsoidDataFromIRMFile("data/cyprus_AMS_groupF.tsv", doNormalize=TRUE)
cyprusMean <- ellMean(cyprusAMS$vector)
sapply(cyprusAMS$logA, ellVolume)
ellVolume(cyprusMean$logA)

# If all of your ellipsoids have the same orientation, then the mean ellipsoid 
# has that same orientation.
synthRot <- rotUniform()
synthLogAs <- replicate(10, sort(rnorm(3)), simplify=FALSE)
synthElls <- lapply(synthLogAs, function(logA) ellEllipsoidFromRotationLogA(synthRot, logA))
synthMean <- ellMean(lapply(synthElls, function(ell) ell$vector))
oriVariance(lapply(synthElls, function(ell) ell$rotation),
            synthMean$rotation, group=oriLineInPlaneGroup)

# Also, if all of your ellipsoids have the same orientation, then the mean 
# ellipsoid's semi-axis lengths are the geometric mean of the ellipsoids' semi-
# axis lengths. Equivalently, the mean ellipsoid's log-semi-axis lengths are 
# the arithmetic mean of the ellipsoids' log-semi-axis lengths.
arithmeticMean(lapply(synthElls, function(ell) ell$logA))
synthMean$logA



### MEASURING DISPERSION ###

# Once the ellipsoids are converted to log-ellipsoid vectors, we can throw all 
# kinds of multivariate statistics at them. For example, we can compute the 
# covariance matrix, and use its eigenvalues to describe the dispersion.

# In the following example, the ellipsoids are un-normalized, so their vectors 
# are 6D, so the covariance matrix is 6x6, so there are six eigenvalues. 
# Notice that most of the eigenvalues are with one order of magnitude, but the 
# greatest eigenvalue is 493 times as big as those.
cyprusAMS <- geoEllipsoidDataFromIRMFile("data/cyprus_AMS_groupF.tsv", doNormalize=FALSE)
ellCovarianceScalars(cyprusAMS$vector)

# After some poking around, we discover that the volumes of the ellipsoids 
# vary by 11 orders of magnitude. Yikes.
hist(sapply(cyprusAMS$logA, ellVolume), xlab="volume")
hist(sapply(cyprusAMS$logA, function(logA) log(ellVolume(logA))), xlab="log(volume)")

# In AMS it is common to volume-normalize the ellipsoids. Doing so reduces 
# them to 5D vectors, whose dispersion is described by five eigenvalues of 
# covariance.
cyprusAMS <- geoEllipsoidDataFromIRMFile("data/cyprus_AMS_groupF.tsv", doNormalize=TRUE)
ellCovarianceScalars(cyprusAMS$vector)

# Covariance is the foundation of principal component analysis, which we'll 
# use seriously in the next tutorial.



### CONCLUSION ###

# When ellipsoids are packaged as log-ellipsoid vectors, lots of statistics, 
# starting with mean and dispersion, are easy.


