


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
# variable in terms of another variable. This tutorial explains scalar-line 
# regression --- that is, fitting a curve of lines parametrized by a scalar. 
# One can use this technique to smooth a data set, predict missing values, etc.



### EXPLORATION ###

# Load some dike data (Scott et al., 2013) from Cyprus. Each dike is tagged 
# with its location (easting and northing, in m) relative to some benchmark 
# location. Plot the stations in map view and the dike poles in equal-area.
ourData <- geoDataFromFile("data/cyprus_inside_corner_field_data.tsv")
plot(x=ourData$easting, y=ourData$northing, xlab="easting (m)", ylab="northing (m)")
lineKambPlot(ourData$pole)

# Also load some dike data from earlier, published maps (various authors). Try 
# mimicking the code above, to plot the data.
mapData <- geoDataFromFile("data/cyprus_inside_corner_map_data.tsv")
# YOUR CODE GOES HERE!

# Assuming that they are compatible with each other (!), combine the two data 
# sets. Plot the dike poles colored by northing, from red (south) to magenta 
# (north). See the vague rainbow pattern? That rainbow suggests that there is 
# a relationship between northing and dike direction.
allData <- rbind(ourData, mapData)
lineEqualAreaPlot(allData$pole, colors=hues(allData$northing))

# Also, in a new window, make a 3D plot showing northing versus equal-area. 
# You can resize, zoom, and rotate the 3D view. This is the same information, 
# visualized in a different way. Try to imagine the best-fit curve through the 
# point cloud.
lineEqualAreaScalarPlot(allData$pole, scales(allData$northing), colors=hues(allData$northing))

# Just so we're clear about what the colors mean, here is the map of stations 
# with the same coloring.
plot(x=allData$easting, y=allData$northing, col=hues(allData$northing), xlab="easting (m)", ylab="northing (m)")



### GEODESIC (GREAT-CIRCLE) REGRESSION ###

# Now we try to fit a simple curve to the dike poles. The curve is a geodesic 
# (great circle arc) on the unit sphere, traversed at a constant rate.
allRegr <- lineRescaledGeodesicRegression(allData$northing, allData$pole, numPoints=100)

# The regression is accomplished using an iterative numerical optimization 
# algorithm. Lots of things can go wrong, so we have to check for errors. 
# First, the $error code should be 0. (If it's not, then we can try increasing 
# the number of iterations allowed.) Second, the $minEigenvalue should be 
# positive. (If it's not, then we can't do much about it.)
allRegr$error
allRegr$minEigenvalue

# Once we're assured that the regression worked, we can start inspecting the 
# results. Here is the curve that it has fit.
lineEqualAreaPlot(allData$pole, colors=hues(allData$northing), curves=list(allRegr$points))

# This is the same plot, but with a curve drawn to connect each datum to its 
# corresponding prediction.
allPreds <- lapply(allData$northing, allRegr$prediction)
allCurves <- thread(lineGeodesicPoints, allData$pole, allPreds)
lineEqualAreaPlot(allData$pole, colors=hues(allData$northing), curves=c(allCurves, list(allRegr$points)))

# Here is the same plot, but emphasizing the predictions rather than the data.
lineEqualAreaPlot(allPreds, colors=hues(allData$northing), curves=c(allCurves, list(allRegr$points)))

# This is the predicted axis of rotation.
geoTrendPlungeDegFromCartesian(allRegr$rotation[,3])

# This is the predicted amount of rotation per m of northing, in degrees. It's 
# a small number. Multiply it by 1000 to get the predicted rotation per km of 
# northing.
allRegr$a / degree
allRegr$a / degree * 1000

# R^2 is a statistic that measures how well the regression describes the data. 
# R^2 is always between 0 and 1. If R^2 = 0, then the regression doesn't 
# explain any of the variation in the data. If R^2 = 1, then the regression 
# perfectly explains the variation in the data (and the data plot on the 
# prediction curve).
allRegr$rSquared

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

# Is our regression result significant, or could it be an artifact of the 
# random variability in the data? To address that question, we perform a 
# permutation test, somewhat like the Wellner (1979) test from the previous 
# tutorial.

# Intuitively, we are wondering whether the relationship between northings and 
# dike pole directions is meaningful. So we arbitrarily reassign the northings 
# to the poles and perform the regression again. If we get a greater R^2, then 
# the new relationship explains the data better than the original one. And the 
# new relationship is meaningless, so maybe the original relationship was 
# meaningless. If we get a lesser R^2, then that's a little evidence that the 
# original relationship had some meaning.

# Concretely, we compute many (say, 1,000 or 10,000) permutations, let p be 
# the proportion in which R^2 exceeds the original R^2, and call the original 
# result significant if p < 0.05, say.

# For the sake of time, we compute only 100 permutations here.
allRSqs <- permutedRSquareds(numPerms=100, lineRescaledGeodesicRegression, allData$northing, allData$pole)

# Because the numerical optimization sometimes fails, you don't always get as 
# many R^2 values as you requested.
length(allRSqs)

# Here is our p-value --- the fraction of successful optimizations that 
# yielded large R^2.
sum(allRSqs > allRegr$rSquared) / length(allRSqs)



### USING THE RESULTS ###

# Besides letting us describe the dike data set in a new way --- 6 degrees of 
# rotation per km of northing --- the regression results give us a couple of 
# insights.

# Recall that our data set has some holes, for example around easting 504000 
# m, northing 3870000 m.
plot(x=allData$easting, y=allData$northing, col=hues(allData$northing), xlab="easting (m)", ylab="northing (m)")

# We can use the regression to predict the dike direction at that location.
geoStrikeDipDegFromCartesian(allRegr$prediction(3870000))

# Also, we now have evidence that the dike poles are not 'independent and 
# identically distributed' (IID), because they depend on geography. Many 
# inference techniques assume that the data are IID. So we should not use 
# those techniques on this data set. That's good to know.



### OTHER IDEAS ###

# Our library offers functions to regress lines along small circles. We just 
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


