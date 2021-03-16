


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

# In statistics, regression is a technique for explaining the variation in one 
# variable in terms of another variable. This tutorial explains scalar-ray 
# regression --- that is, fitting a curve of rays parametrized by a scalar. 
# One can use this technique to smooth a data set, predict missing values, etc.



### EXPLORATION ###

# Let's again load a data set of paleomagnetic directions from Cyprus, 
# compiled from various papers (Bonhommet et al., 1988; Allerton, 1989; 
# Abelson et al., 2002; Granot et al., 2006; Scott et al., 2013).
pmag <- geoDataFromFile("data/cyprusSouthInPmag.csv")
rayEqualAreaPlot(pmag$direction)

# Each of the data has a geographic location attached to it, expressed as 
# easting and northing (in m) relative to a benchmark. Here are the stations 
# in map view.
plot(x=pmag$easting, y=pmag$northing, xlab="easting (m)", ylab="northing (m)")

# Playing around with the data, we plot them colored by easting. Notice the 
# vague rainbow in this plot: from reds, oranges, and yellows, to greens, 
# blues, and magentas. This rainbow suggest that the paleomagnetic directions 
# depend on easting.
rayEqualAreaPlot(pmag$direction, colors=hues(pmag$easting))

# Just so we're clear about what the colors mean, here are the stations 
# colored in the same way.
plot(x=pmag$easting, y=pmag$northing, col=hues(pmag$easting), xlab="easting (m)", ylab="northing (m)")

# (Actually, Sarah suspects that the paleomagnetic directions depend on 
# distance from a certain fossil spreading ridge. The ridge strikes 010 
# degrees. For the sake of simplicity in this tutorial, we call the ridge NS, 
# so that distance from the ridge is essentially easting.)



### GEODESIC (GREAT-CIRCLE) REGRESSION ###

# Now we try to fit a simple curve to the paleomagnetic data. The curve is a 
# geodesic (great circle arc) on the unit sphere, traversed at a constant rate.
pmagRegr <- rayRescaledGeodesicRegression(pmag$easting, pmag$direction, numPoints=100)

# The regression is accomplished using an iterative numerical optimization 
# algorithm. Lots of things can go wrong, so we have to check for errors. 
# First, the $error code should be 0. (If it's not, then we can try increasing 
# the number of iterations allowed.) Second, the $minEigenvalue should be 
# positive. (If it's not, then we can't do much about it.)
pmagRegr$error
pmagRegr$minEigenvalue

# Once we're assured that the regression worked, we can start inspecting the 
# results. Here is the curve that it has fit.
rayEqualAreaPlot(pmag$direction, colors=hues(pmag$easting), curves=list(pmagRegr$points))

# This is the same plot, but with a curve drawn to connect each datum to its 
# corresponding prediction.
pmagPreds <- lapply(pmag$easting, pmagRegr$prediction)
pmagCurves <- thread(rayGeodesicPoints, pmag$direction, pmagPreds)
rayEqualAreaPlot(pmag$direction, colors=hues(pmag$easting), curves=c(pmagCurves, list(pmagRegr$points)))

# Here is the same plot, but emphasizing the predictions rather than the data.
rayEqualAreaPlot(pmagPreds, colors=hues(pmag$easting), curves=c(pmagCurves, list(pmagRegr$points)))

# This is the predicted axis of rotation.
geoTrendPlungeDegFromCartesian(pmagRegr$rotation[,3])

# This is the predicted amount of rotation per m of easting, in degrees. It's 
# a small number. Multiply it by 1000 to get the predicted rotation per km of 
# easting.
pmagRegr$a / degree
pmagRegr$a / degree * 1000

# R^2 is a statistic that measures how well the regression describes the data. 
# R^2 is always between 0 and 1. If R^2 = 0, then the regression doesn't 
# explain any of the variation in the data. If R^2 = 1, then the regression 
# perfectly explains the variation in the data (and the data plot on the 
# prediction curve).
pmagRegr$rSquared

# Which values of R^2 qualify as 'good'? The answer depends on your 
# expectations. If you're a particle physicist, then you expect your data to 
# be explainable by simple and precise laws, so you might not be satisfied 
# with R^2 = 0.8. If you're a stock market analyst, then you expect the 
# motions of stock prices to be unexplainable, so you might be thrilled to 
# discover a tendency with R^2 = 0.1.

# In our experience with geologic data sets, R^2 values between 0.2 and 0.7 
# are pretty typical. Rocks are governed by physical laws, but they are also 
# complicated and messy, and there is no reason to believe that a simple 
# regression can describe their behavior precisely.



### SIGNIFICANCE ###

# A separate issue is whether the regression results are statistically 
# significant. In other words, could the tendency that we've discovered be the 
# result of random variation in the data, rather than a meaningful geographic 
# signal?

# To address this issue, we perform a permutation test. We randomly permute 
# the eastings (but not the paleomagnetic directions), re-run the regression, 
# and compute its R^2. We do this many times, to get a bunch of R^2 values. 
# For the sake of time, we'll do just 100 permutations here. For a real 
# research problem you would want to do at least 1,000 permutations.
pmagRSqs <- permutedRSquareds(numPerms=100, rayRescaledGeodesicRegression, pmag$easting, pmag$direction)

# Because the numerical optimization sometimes fails, you don't always get as 
# many R^2 values as you requested.
length(pmagRSqs)

# Now here's the interesting part. What if a random permutation of the 
# eastings yielded an R^2 greater than the original R^2? That would be 
# evidence that the original relationship between easting and paleomagnetic 
# direction was meaningless. So let p be the fraction of permutations that 
# yield larger R^2 values than the original. Small values of p (say, p < 0.05) 
# indicate that the original relationship between easting and direction is 
# meaningful, and hence that the result is significant.
sum(pmagRSqs > pmagRegr$rSquared) / length(pmagRSqs)



### USING THE RESULTS ###

# Besides letting us describe the paleomagnetic data set in a new way --- 
# about 2 degrees of rotation per km of easting --- the regression results 
# give us a couple of insights.

# Recall that our data set has some holes, for example around easting 505000 
# m, northing 3860000 m.
plot(x=pmag$easting, y=pmag$northing, col=hues(pmag$easting), xlab="easting (m)", ylab="northing (m)")

# We can use the regression to predict the paleomagnetic direction at that 
# location.
geoTrendPlungeDegFromCartesian(pmagRegr$prediction(505000))

# Also, we now have evidence that the paleomagnetic directions are not 
# 'independent and identically distributed' (IID), because they depend on 
# geography. Many inference techniques assume that the data are IID. So we 
# should not use those techniques on this data set. That's good to know.



### OTHER IDEAS ###

# Our library offers functions to regress rays along small circles. We just 
# haven't written tutorials for them yet.

# Another technique we've used is kernel regression. You can see it in our 
# rotation tutorials.

# The problem with all of these purely statistical techniques is that they are 
# not motivated by any underlying geologic theory. An approach called 'inverse 
# modeling' can be used to fit geologically informed models. It's beyond the 
# scope of these tutorials.



### CONCLUSION ###

# In geology, it is common for data to exhibit spatial or temporal patterns. 
# Regression is one tool for quantifying such patterns. But you need to make 
# sure that you're not over-detecting patterns that aren't really there.


