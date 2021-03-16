


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
# dispersed or 'spread-out' about that mean the data are. This tutorial 
# illustrates the analogous statistics for rays. It turns out that there are 
# several options for quantifying mean and dispersion.



### EXTRINSIC METHODS ###

# Here is a set of paleomagnetic directions from Cyprus, compiled from various 
# papers (Bonhommet et al., 1988; Allerton, 1989; Abelson et al., 2002; Granot 
# et al., 2006; Scott et al., 2013).
pmag <- geoDataFromFile("data/cyprusSouthInPmag.csv")
rayEqualAreaPlot(pmag$direction)

# Here's the simplest idea for averaging a bunch of rays r1, r2, ..., rn: Just 
# take their arithmetic mean arith := (r1 + r2 + ... + rn) / n. The problem is 
# that rays are encoded mathematically as unit vectors (or, equivalently, 
# points on the unit sphere), and arith is usually not a unit vector. But we 
# can divide it by its length to project it back to the unit sphere. This 
# concept is therefore called the 'projected arithmetic mean' or 'extrinsic 
# mean'.
pmagMean <- rayProjectedMean(pmag$direction)
rayEqualAreaPlotTwo(pmag$direction, list(pmagMean), colorB="red")
pmagMean
geoTrendPlungeDegFromCartesian(pmagMean)

# The quantity 1 - |arith| measures the dispersion of the data set on a scale 
# from 0 to 1. You can get it through this function.
pmagMeanScat <- rayMeanScatter(pmag$direction)
pmagMeanScat$scatter

# Here are two important end members. First, when all of the rays ri are 
# identical, then they also equal arith, which thus has length 1, so that the 
# dispersion is 0. This highly concentrated synthetic example is close to this 
# end member.
synth <- rayFisher(rayUniform(), kappa=500, 100)
rayEqualAreaPlot(synth)
rayMeanScatter(synth)$scatter

# Second, when all of the rays ri point in 'perfectly different' directions, 
# arith is the zero vector, so the dispersion is 1 (and the projected 
# arithmetic mean is undefined). This nearly uniform synthetic example is 
# close to this end member.
synth <- rayUniform(1000)
rayEqualAreaPlot(synth)
rayMeanScatter(synth)$scatter



### DISTRIBUTIONAL METHODS ###

# In elementary statistics, the normal distribution plays a pivotal role. In 
# fact, computing the sample mean and standard deviation is tantamount to 
# fitting a normal distribution to the data. In statistics of rays, there are 
# at least two versions of the normal distribution, called Fisher and Kent.

# The Fisher distribution is simple, in that it is isotropic (rotationally 
# symmetric) about its mean. If we assume that the data arise from a Fisher 
# distribution, then this function estimates the mean (muHat) and 
# concentration (kappaHat) of that distribution. kappaHat near 0 means the 
# data are very dispersed; kappaHat near infinity means the data are very 
# concentrated. Our paleomagnetic data have kappaHat near 4.6.
rayFisherMLE(pmag$direction)

# It's noteworthy that the mean produced by this technique is always identical 
# to the projected arithmetic mean above.

# Unfortunately, the data don't look very isotropic about their mean, so the 
# Fisher distribution is a poor choice. The Kent distribution is capable of 
# richer behavior including anisotropy about the mean. There are some Kent 
# tools in the R package Directional. They have not yet been integrated into 
# our geologyGeometry library.



### INTRINSIC METHODS ###

# The variance of the data {r1, r2, ..., rn} about a ray r is proportional to 
# the sum of the squared distances from r to the ri. The following code 
# defines a function to help you visualize this concept. You can then give 
# that function any trend and plunge you want, and observe how the variance 
# changes.
varianceVisualizer <- function(r) {
  rayEqualAreaPlot(pmag$direction, curves=lapply(pmag$direction, rayGeodesicPoints, r))
  rayVariance(pmag$direction, r)
}
varianceVisualizer(geoCartesianFromTrendPlungeDeg(c(220, 45)))

# The Frechet mean is defined to be the ray r that minimizes the variance.
pmagMeanVar <- lineMeanVariance(pmag$direction)
varianceVisualizer(pmagMeanVar$mean)
geoTrendPlungeDegFromCartesian(pmagMeanVar$mean)

# The projected arithmetic mean (green in plot below) is usually close to, but 
# not identical to, the Frechet mean (blue in plot).
rayEqualAreaPlotThree(pmag$direction, list(pmagMean), list(pmagMeanVar$mean), colorB="green", colorC="blue")
rayDistance(pmagMean, pmagMeanVar$mean) / degree

# The variance of the data about their Frechet mean is a scalar measure of 
# dispersion.
pmagMeanVar$variance

# Another approach is to transfer the data into the tangent plane to the unit 
# sphere at the mean, and do principal component analysis in that plane. For 
# starters, PCA gives two 'magnitudes' that characterize the dispersion in the 
# data set anisotropically.
pmagPCA <- rayPrincipalComponentAnalysis(pmag$direction, pmagMeanVar$mean, 10)
pmagPCA$magnitudes

# Corresponding to the two magnitudes are two directions. The first direction 
# captures girdling, and the second is perpendicular to the first.
rayEqualAreaPlot(pmag$direction, curves=pmagPCA$curves)

# You can imagine those two curves as x- and y-axes in the tangent plane. Here 
# are the data in the tangent plane, with that choice of axes. Many 
# statistical techniques can be applied to the data in this format (but let's 
# not).
xys <- sapply(pmag$direction, pmagPCA$pcsFromRay)
plot(x=xys[1,], y=xys[2,])

# When we transfer a data set from a sphere to a tangent plane, we inevitably 
# distort the geometry of that data set --- especially for data far from the 
# point of tangency. So this PCA technique should be used only for tightly 
# concentrated data sets (more concentrated than this example, actually).



### WHEN NOT TO COMPUTE THE MEAN ###

# Here is a paleomagnetic data set of Chapman et al. (in prep) from Iceland. 
# Notice that there seems to be some polarity reversal in the data.
icePmag <- geoDataFromFile("data/icelandOurFlowInSitu29March.csv")
rayEqualAreaPlot(icePmag$direction)

# In the following plot, the projected arithmetic mean is shown in green, and 
# the Frechet mean is shown in blue. Notice how far apart they are. Also 
# notice that neither mean is really 'in' the point cloud.
iceProj <- rayProjectedMean(icePmag$direction)
iceFrechet <- rayMeanVariance(icePmag$direction)$mean
rayEqualAreaPlotThree(icePmag$direction, list(iceProj), list(iceFrechet), colorA="black", colorB="green", colorC="blue")

# The means are behaving badly because the data set is bimodal: There are two 
# regions of high density. Each mean is trying to describe both regions but 
# instead falling between them and describing neither. The mean is just not a 
# good summary statistic for a bimodal data set. (And measures of dispersion 
# about the mean are also a bad idea).

# In this example, one option is to forget polarity and treat the data as mere 
# lines. Another option is to isolate the two modes and analyze them 
# separately. That task can be accomplished using clustering algorithms or 
# mixtures of distributions. In this tutorial, let's simply isolate the 
# downward- and upward-pointing rays.

# Here are the downward-pointing rays and their means. The means are 
# behaving much better than they did above.
icePmagDown <- Filter(function(u) {u[[3]] <= 0}, icePmag$direction)
iceProjDown <- rayProjectedMean(icePmagDown)
iceFrechetDown <- rayMeanVariance(icePmagDown)$mean
rayEqualAreaPlotThree(icePmagDown, list(iceProjDown), list(iceFrechetDown), colorA="black", colorB="green", colorC="blue")

# Similarly, here are the upward-pointing rays and their means.
icePmagUp <- Filter(function(u) {u[[3]] > 0}, icePmag$direction)
iceProjUp <- rayProjectedMean(icePmagUp)
iceFrechetUp <- rayMeanVariance(icePmagUp)$mean
rayEqualAreaPlotThree(icePmagUp, list(iceProjUp), list(iceFrechetUp), colorA="black", colorB="green", colorC="blue")



### CONCLUSION ###

# To describe your data set, you typically want to plot the data, compute the 
# mean, and compute some measure of dispersion. You have various options.


