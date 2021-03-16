


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

# Given two ellipsoidal data sets about the same rocks, it is natural to ask 
# whether they tell us the same thing. In this tutorial we do an example: a 
# paired two-sample hypothesis test, which is really a kind of one-sample 
# hypothesis test.



### EXPLORATION ###

# Load a data set from New Caledonia (Titus et al., 2011; Chatzaras et al., in 
# prep). The data have been pared down to just 11 field sites where we have 
# two measures of shape preferred orientation (SPO). The first SPO comes from 
# field measurements of orthopyroxene, with the SPO ellipsoid reconstructed 
# according to the method of Robin (2002). The second SPO comes from averaging 
# spinel grains measured by X-ray computed tomography.
source("data/newcalOPXSpinelSPO.R")
length(ncOrthoSpinels)

# Inspect the first datum, just to see what it's like. It's a station name, an 
# easting-northing in km relative to the supposed ridge-transform 
# intersection, an orthopyroxene ellipsoid, and a spinel ellipsoid.
ncOrthoSpinels[[1]]

# Here are the orthopyroxene shapes (red) and spinel shapes (blue), drawn with 
# lines connecting each pair. The orthopyroxene ellipsoids seem closer to 
# spherical than the spinel ellipsoids are.
ncOrthoLogAs <- lapply(ncOrthoSpinels, function(station) station$ortho$logA)
ncSpinelLogAs <- lapply(ncOrthoSpinels, function(station) station$spinel$logA)
ncLogAPairs <- lapply(1:length(ncOrthoSpinels), function(i) list(ncOrthoLogAs[[i]], ncSpinelLogAs[[i]]))
ellHsuNadaiPlot(unlist(ncLogAPairs, recursive=FALSE, use.names=FALSE),
                curves=ncLogAPairs, colors=c("red", "blue"))

# Here is a similar plot of the orientations. I don't see much of a pattern.
ncOrthoRotations <- lapply(ncOrthoSpinels, function(station) station$ortho$rotation)
ncSpinelRotations <- lapply(ncOrthoSpinels, function(station) station$spinel$rotation)
ncRotationPairs <- lapply(1:length(ncOrthoSpinels), function(i) list(ncOrthoRotations[[i]], ncSpinelRotations[[i]]))
ellEqualVolumePlot(unlist(ncRotationPairs, recursive=FALSE, use.names=FALSE),
                   unlist(ncLogAPairs, recursive=FALSE, use.names=FALSE),
                   rotCurves=ncRotationPairs, aCurves=ncLogAPairs, colors=c("red", "blue"))

# And here's a similar plot of the ellipsoid vectors. The orthopyroxene 
# ellipsoids are closer to the origin and hence closer to spherical.
ncOrthoVectors <- lapply(ncOrthoSpinels, function(station) station$ortho$vector)
ncSpinelVectors <- lapply(ncOrthoSpinels, function(station) station$spinel$vector)
ncVectorPairs <- lapply(1:length(ncOrthoSpinels), function(i) list(ncOrthoVectors[[i]], ncSpinelVectors[[i]]))
ellVectorPlot(c(1, 2, 3), unlist(ncVectorPairs, recursive=FALSE, use.names=FALSE),
              curves=ncVectorPairs, colors=c("red", "blue"))

# Here are the ellipsoid pairwise plots. Right now this kind of plot doesn't 
# let us draw the connecting lines. But it's again clear that the 
# orthopyroxene ellipsoids are generally closer to the origin and hence more 
# spherical.
ellPairsPlot(unlist(ncVectorPairs, recursive=FALSE, use.names=FALSE), colors=c("red", "blue"))



### PAIRED HYPOTHESIS TEST ###

# We're going to do a two-sample hypothesis test --- that is, a test involving 
# two data sets. But this two-sample test is special, because the data are 
# paired: Each orthopyroxene ellipsoid has a corresponding spinel ellipsoid. 
# We're wondering whether they're the same, which is the same as wondering 
# whether their differences are zero. Look closely at this plot of 
# differences. The point cloud does not obviously miss the origin.
ncDiffVectors <- lapply(1:length(ncOrthoSpinels), function(i) {ncSpinelVectors[[i]] - ncOrthoVectors[[i]]})
ellPairsPlot(ncDiffVectors)

# A paired two-sample hypothesis test is very easy. It's just a one-sample 
# hypothesis test on the differences. So, mimicking our earlier ellipsoid 
# inference tutorial, we do this.
ncInfT2 <- ellHotellingT2Inference(ncDiffVectors, c(0, 0, 0, 0, 0))
ncInfMM <- ellBootstrapMMInference(ncDiffVectors, c(0, 0, 0, 0, 0), numBoots=1000)
ncInfS <- ellBootstrapSInference(ncDiffVectors, c(0, 0, 0, 0, 0), numBoots=1000)
ncInf <- ellBootstrapInference(ncDiffVectors, numBoots=1000)

# Sometimes I get an error on the MM-type test or the S-type test. I think 
# that it's a bug in those routines. If that happens to you, then ignore the 
# corresponding lines below. Here are the p-values.
ncInfT2$p.value
ncInfMM$p.value
ncInfS$p.value
ncInf$pvalue(c(0, 0, 0, 0, 0))

# The T2 test does not reject the null hypothesis. I'm not really sure that 
# the sample size was big enough to justify use of that method anyway. I don't 
# trust the MM- and S-type tests here. I've tried the last bootstrap test a 
# few times, and I always get p-values near 0.05. So let's crank up the number 
# of bootstraps.
ncInf <- ellBootstrapInference(ncDiffVectors, numBoots=100000)
ncInf$pvalue(c(0, 0, 0, 0, 0))

# With 100,000 bootstraps, I tend to get p-values of 0.050 or 0.051. So I 
# can't reject the null hypothesis. The test is inconclusive. We do not know 
# that field orthopyroxene SPO and XRCT spinel SPO are different.

# Nor do we know that they're the same. But if we had some other reason to 
# believe that they're the same, then it would be reasonable to report the 
# hypothesis test results and proceed under that assumption.



### OTHER IDEAS ###

# The preceding analysis implicitly assumes that the differences are 
# independent and identically distributed. But maybe they're not. Maybe 
# something about the field location, orthopyroxene SPO, or spinel SPO 
# correlates with the difference. We could explore that possibility using 
# regression, for example.



### CONCLUSION ###

# A two-sample test can help us understand whether two populations are 
# different. In a paired two-sample test, the data from one population 
# correspond one-to-one with the data from the other population, and we test 
# their difference.


