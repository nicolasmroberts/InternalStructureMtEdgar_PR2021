


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

# In rays, lines, and rotations, there are three broad approaches to 
# quantifying the mean and dispersion of a data set: extrinsic (projected), 
# distributional (e.g., Fisher, Kent, Watson, Bingham, matrix Fisher), and 
# intrinsic (Frechet). In orientations, the symmetry group complicates the 
# calculations, so that the extrinsic approach doesn't work and distributions 
# are a bit finicky. So our best option is the intrinsic methods: Frechet 
# mean, variance, tangent space approximation, etc.



### FRECHET MEAN AND VARIANCE ###

# Load a set foliation-lineation orientations from the western Idaho shear 
# zone (Giorgis and Tikoff, 2004). Because of the 4-fold symmetry, the data 
# appear four times in a rotation plot.
wiszData <- geoDataFromFile("data/wiszFollins.tsv")
oriEqualAnglePlot(wiszData$rotation, group=oriLineInPlaneGroup)

# The distance between two orientations is defined to be the size of the 
# smallest rotation that rotates one to the other. Beyond that, the variance 
# and the Frechet mean are defined exactly as they are for rays, lines, and 
# rotations. In the following plot, the mean is at the center of the spider, 
# and the variance is proportional to the sum of the squared lengths of the 
# spider legs.
wiszMeanVar <- oriMeanVariance(wiszData$rotation, group=oriLineInPlaneGroup)
wiszCurves <- lapply(
  oriNearestRepresentatives(wiszData$rotation, wiszMeanVar$mean, group=oriLineInPlaneGroup),
  rotGeodesicPoints, wiszMeanVar$mean)
oriEqualAnglePlot(wiszData$rotation, curves=wiszCurves, group=oriLineInPlaneGroup)

# The variance of a data set about its Frechet mean is a scalar measure of the 
# dispersion.
wiszMeanVar$variance



### COMPARED TO LINE MEANS SEPARATELY ###

# Remember that analyzing the foliations and lineations separately is 'wrong'
# (or not ideal, at least) because it ignores information present in the data 
# set --- namely, how they are paired and correlated. But you might be curious 
# how wrong it is.

# Here are the mean foliation pole and the mean lineation direction. They're 
# not perpendicular, but they're very close to it.
wiszMeanFol <- lineProjectedMean(wiszData$pole)
wiszMeanLin <- lineProjectedMean(wiszData$direction)
lineDistance(wiszMeanFol, wiszMeanLin) / degree

# In this plot, the foliation pole arising from the orientation mean is in 
# green, and the mean of the foliation poles is in blue. They're close.
lineEqualAreaPlotThree(wiszData$pole, list(wiszMeanVar$mean[1,]), list(wiszMeanFol), colorB="green", colorC="blue")

# Similarly, the lineation arising from the orientation mean is close to the 
# mean of the lineations.
lineEqualAreaPlotThree(wiszData$direction, list(wiszMeanVar$mean[2,]), list(wiszMeanLin), colorB="green", colorC="blue")

# So, there doesn't seem to be much harm in doing things the wrong way. But 
# consider: In other examples the harm might be greater. And in a more 
# complicated analysis, even small discrepancies, like the ones observed in 
# this example, can accumulate to produce big problems. So why not just do it 
# the right way, if you can?



### WHEN NOT TO COMPUTE THE MEAN ###

# In rays, lines, and rotations we've seen that multi-modal data sets are not 
# well summarized by the mean. The same is true for orientations. However, the 
# symmetry group makes detecting multiple modes in a plot a little trickier.

# Here's a synthetic data set. Is it multi-modal?
synthData <- rotIsotropicFisher(rotMatrixFromAxisAngle(c(1, 0, 0, pi)), kappa=10, n=100)
oriEqualVolumePlot(synthData, group=oriTrivialGroup)

# Here's a synthetic data set. Is it multi-modal?
synthData <- c(
  rotIsotropicFisher(rotMatrixFromAxisAngle(c(1, 0, 0, pi)), kappa=10, n=100),
  rotIsotropicFisher(rotMatrixFromAxisAngle(c(0, 1, 0, pi)), kappa=10, n=100))
oriEqualVolumePlot(synthData, group=oriLineInPlaneGroup)

# Here's a synthetic data set. Is it multi-modal?
synthData <- c(
  rotIsotropicFisher(rotMatrixFromAxisAngle(c(1, 0, 0, pi)), kappa=10, n=100),
  rotIsotropicFisher(rotMatrixFromAxisAngle(c(0, 1, 0, pi)), kappa=10, n=100))
oriEqualVolumePlot(synthData, group=oriRayInPlaneGroup)



### CONCLUSION ###

# Because of the symmetry group, techniques for orientations are more limited 
# than techniques for rays, lines, and rotations. In particular, I know of 
# only one good way to conceptualize the mean.


