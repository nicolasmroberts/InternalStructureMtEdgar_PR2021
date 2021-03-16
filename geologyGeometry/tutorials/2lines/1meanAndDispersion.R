


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

# This tutorial illustrates the analogous statistics for lines. It turns out 
# that there are several options for quantifying mean and dispersion.



### EXTRINSIC METHODS ###

# Let's make sure our dike poles from Cyprus site 230 (Titus, pers. comm.) are 
# loaded.
site230 <- geoDataFromFile("data/cyprusDikesSite230.tsv")
lineEqualAreaPlot(site230$pole)

# Here is one notion of mean. It is constructed by forming a 'scatter matrix' 
# and then computing the eigenvector with greatest eigenvalue.
site230Mean <- lineProjectedMean(site230$pole)
lineEqualAreaPlotTwo(site230$pole, list(site230Mean), colorB="red")

# Here's the mean in various numeric formats.
site230Mean
geoTrendPlungeDegFromCartesian(site230Mean)
geoStrikeDipDegFromCartesian(site230Mean)

# The scatter matrix contains more information than just the mean. To get at 
# this information, we use a different function: lineMeanScatter. The 
# eigenvectors of the scatter matrix are three perpendicular lines. The first 
# line is the mean. The second line defines the direction of greatest girdling 
# from the mean. The third line is the pole to the girdle.
site230Scatter <- lineMeanScatter(site230$pole)
lineEqualAreaPlotTwo(site230$pole,
                     list(site230Scatter$vectors[,1], site230Scatter$vectors[,2], site230Scatter$vectors[,3]),
                     colorB="red", curves=list(rayGreatCircle(site230Scatter$vectors[,3])))

# The eigenvalues of the scatter matrix are non-negative and sum to 1. 
# Roughly, the first number quantifies how concentrated the data are (about 
# the mean). The second number quantifies how girdled the data are (away from 
# the mean). The third number quantifies how spread-out-over-the-sphere the 
# data are (away from the girdle).
site230Scatter$values



### DISTRIBUTIONAL METHODS ###

# In elementary statistics, the normal distribution plays a pivotal role. In 
# fact, computing the sample mean and standard deviation is tantamount to 
# fitting a normal distribution to the data. In statistics of lines, there are 
# at least two versions of the normal distribution, called Watson and Bingham.

# The Watson distribution is simple, in that it is isotropic (rotationally 
# symmetric) about its mean. If we assume that the data arise from a Watson 
# distribution, then this function estimates the mean (muHat) and 
# concentration (kappaHat) of that distribution. The mean is identical to the 
# scatter-matrix mean above. The inverse of kappaHat is a measure of 
# dispersion.
lineWatsonMLE(site230$pole)

# Unfortunately, the data don't look very isotropic about their mean, so the 
# Watson distribution might be a poor choice. Let's try Bingham. Our library 
# offers two methods. The first one is reliable but slow and a little 
# hands-on. First you check that 'error' is 0 and 'minEigenvalue' is positive. 
# (If not, then the calculation went badly.)
site230MLE <- lineBinghamMLE(site230$pole)
site230MLE$error
site230MLE$minEigenvalue

# Then the 'vectors' describe the directions of dispersion. Actually they're 
# identical to the scatter-matrix eigenvectors above. In particular, the first 
# eigenvector is the scatter-matrix mean.
site230MLE$vectors[,1]

# The 'values' always add to 0. Their negatives describe the amount of 
# concentration.
site230MLE$values

# Thus the inverses of the negatives of the eigenvalues describe the amount of 
# dispersion. Admittedly this is all a bit cryptic.
1 / (-site230MLE$values)

# A lesson here is: Because the Bingham distribution is richer than the Watson 
# distribution, it is able to describe data more subtly. However, the 
# description is so subtle that it takes a lot of training to decipher it.

# Anyway, here is the second method provided by our library. Because it uses 
# some approximations, it's faster than the one above. But sometimes it just 
# doesn't give an answer --- in fact, this time.
lineTauxeBinghamMLE(site230$pole)



### INTRINSIC METHODS ###

# The variance of the data {l1, l2, ..., ln} about a line l is proportional to 
# the sum of the squared distances from l to the li. The following code 
# defines a function to help you visualize this concept. You can then give 
# that function any trend and plunge you want, and observe how the variance 
# changes.
varianceVisualizer <- function(l) {
  lineEqualAreaPlot(site230$pole, curves=lapply(site230$pole, lineGeodesicPoints, l))
  lineVariance(site230$pole, l)
}
varianceVisualizer(geoCartesianFromTrendPlungeDeg(c(220, 45)))

# The Frechet mean is defined to be the line l that minimizes the variance.
site230MeanVar <- lineMeanVariance(site230$pole)
varianceVisualizer(site230MeanVar$mean)
geoTrendPlungeDegFromCartesian(site230MeanVar$mean)

# The projected arithmetic mean (green in plot below) is usually close to, but 
# not identical to, the Frechet mean (blue in plot).
lineEqualAreaPlotThree(site230$pole, list(site230Mean), list(site230MeanVar$mean), colorB="green", colorC="blue")
lineDistance(site230Mean, site230MeanVar$mean) / degree

# The variance of the data about their Frechet mean is a scalar measure of 
# dispersion.
site230MeanVar$variance



### WHEN NOT TO COMPUTE THE MEAN ###

# The mean is a good summary of the 'center' of a data set only when the data 
# set has a meaningful center. To illustrate this point, let's randomly 
# generate a synthetic data set by sampling from two Watson distributions.
synthFirst <- lineWatson(mu=lineUniform(), kappa=30, n=50)
synthSecond <- lineWatson(mu=lineUniform(), kappa=30, n=100)
synth <- c(synthFirst, synthSecond)
lineEqualAreaPlot(synth)

# Probably your equal-area plot shows two distinct regions of high density. 
# (If not, then re-generate the data set until it does.) We say that the data 
# are 'bimodal'. What happens if we compute the mean? The mean falls between 
# the two modes, summarizing neither of them well.
synthProj <- lineProjectedMean(synth)
lineEqualAreaPlotTwo(synth, list(synthProj), colorB="green")

# Incidentally, the projected arithmetic mean and the Frechet mean often 
# differ greatly for bimodal data sets.
synthFrechet <- lineMeanVariance(synth)$mean
lineEqualAreaPlotThree(synth, list(synthProj), list(synthFrechet), colorB="green", colorC="blue")
lineDistance(synthProj, synthFrechet) / degree

# Perhaps we should split the data set into its two modes and analyze them 
# separately. Clustering and mixed distributions are two techniques for 
# performing this splitting, outside the scope of this tutorial.



### CONCLUSION ###

# To describe your data set, you typically want to plot the data, compute the 
# mean, and compute some measure of dispersion. You have various options.


