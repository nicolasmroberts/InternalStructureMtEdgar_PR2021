


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
# variable in terms of another variable. This tutorial explains scalar-
# orientation regression --- that is, fitting a curve of orientations 
# parametrized by a scalar. One can use this technique to smooth a data set, 
# predict missing values, etc.



### EXPLORATION ###

# Load some foliation-lineation orientations from the Ahsahka shear zone, 
# Idaho (Stetson-Lee, 2015; Roberts et al., in review). They are tagged with 
# their geographic locations, expressed as eastings and northings (in m) 
# relative to some benchmark location. Plot them in map view.
ahsahka <- geoDataFromFile("data/Follins_Ahsahka.csv")
plot(x=ahsahka$easting, y=ahsahka$northing, xlab="easting (m)", ylab="northing (m)")

# The Ahsahka shear zone strikes approximately NW-SE. So northeasting is a 
# distance measured perpendicular to the shear zone. We plot the orientations 
# colored by northeasting. A pattern emerges: The reds and yellows plot 
# together, and the blues and purples plot together.
ahsahka$ne <- (ahsahka$easting + ahsahka$northing) / sqrt(2)
oriEqualVolumePlot(ahsahka$rotation, group=oriLineInPlaneGroup, colors=hues(ahsahka$ne))



### GEODESIC REGRESSION ###

# We try to quantify this pattern using a simple kind of regression. This 
# regression fits a curve to the orientation data, based on northeasting, such 
# that there is a steady amount of rotation per m of northeasting.
ahsahkaRegr <- oriRescaledGeodesicRegression(ahsahka$ne, ahsahka$rotation, group=oriLineInPlaneGroup, numSteps=1000)

# The regression relies on an iterative numerical optimization procedure. Lots 
# of things can go wrong. We have to check that the 'error' equals 0. (If not, 
# then try increasing numSteps above.) We also must check that the 
# 'minEigenvalue' is positive. (If not, there's not much we can do.)
ahsahkaRegr$error
ahsahkaRegr$minEigenvalue

# For each value of northeasting, the regression predicts a particular 
# orientation. In other words, it predicts a curve's worth of orientations. 
# You can see that the curve travels from the reds to the blues, approximately.
ahsahkaCurve <- lapply(seq(from=min(ahsahka$ne), to=max(ahsahka$ne), length.out=100), ahsahkaRegr$prediction)
oriEqualVolumePlot(ahsahka$rotation, group=oriLineInPlaneGroup, colors=hues(ahsahka$ne), curves=list(ahsahkaCurve), curveWidth=3)

# Here's the same plot, but with additional curves connecting each datum to 
# where the curve predicts it should be.
ahsahkaPreds <- lapply(ahsahka$ne, ahsahkaRegr$prediction)
ahsahkaCurves <- thread(function(pred, obs) rotGeodesicPoints(pred, oriNearestRepresentative(obs, pred, group=oriLineInPlaneGroup)),
                        ahsahkaPreds, ahsahka$rotation)
oriEqualVolumePlot(ahsahka$rotation, group=oriLineInPlaneGroup, colors=hues(ahsahka$ne), 
                   curves=c(ahsahkaCurves, list(ahsahkaCurve)), curveWidth=3)

# Here's the same plot, but emphasizing the predictions instead of the data.
oriEqualVolumePlot(ahsahkaPreds, group=oriLineInPlaneGroup, colors=hues(ahsahka$ne), 
                   curves=c(ahsahkaCurves, list(ahsahkaCurve)), curveWidth=3)

# The R^2 statistic measures the explanatory power of the regression. It 
# varies between 0 and 1. A value of R^2 = 0 means that the regression 
# explains none of the variation in the data. A value of R^2 = 1 means that 
# the regression explains all of the variation (and the data lie perfectly on 
# the regression curve).
ahsahkaRegr$rSquared

# Which values of R^2 are good enough varies from discipline to discipline. In 
# our experience with geologic data sets, R^2 between 0.2 and 0.7 are typical.



### SIGNIFICANCE ###

# The tendency that we just discovered is not very strong. But a separate 
# issue is: Is it meaningful? Have we discovered a real signal in the data, or 
# could this result have arisen from random variation that has nothing to do 
# with northeasting?

# To address this question, we perform a permutation test. In this test, we 
# randomly permute the northeastings (but not the orientations), re-run the 
# regression, and compute its R^2. Because the northeastings are randomized, 
# they have no meaningful relationship to the orientations. So if this new 
# regression produces a larger R^2 than the original R^2, then we can take 
# that as evidence that the original relationship between northeasting and 
# orientation was also not meaningful.

# We repeat this process many (say, 1,000 or 10,000) times, to get a big list 
# of R^2 values. Let p be the fraction of those values that are greater than 
# the original R^2. Small values of p, such as p < 0.05, indicate a 
# statistically significant result.

# In this tutorial, for the sake of time, we perform only 100 permutations. 
# And even that takes several minutes. If you get tired of waiting, then press 
# RStudio's 'Stop' button to cancel the computation.
ahsahkaRSqs <- permutedRSquareds(numPerms=100, oriRescaledGeodesicRegression, ahsahka$ne, 
                                 ahsahka$rotation, group=oriLineInPlaneGroup, numSteps=10000)

# Because the numerical optimization sometimes fails, you don't always get as 
# many R^2 values as you requested.
length(ahsahkaRSqs)

# Here is our p-value --- the fraction of successful optimizations that 
# yielded large R^2.
sum(ahsahkaRSqs > ahsahkaRegr$rSquared) / length(ahsahkaRSqs)



### USING THE RESULTS ###

# Here's a plot of foliation poles (circles) and lineation directions 
# (squares), colored by northeasting, with the regression prediction curves 
# superimposed.
ahsahkaFols <- lapply(ahsahkaCurve, function(r) r[1,])
ahsahkaLins <- lapply(ahsahkaCurve, function(r) r[2,])
lineEqualAreaPlot(c(ahsahka$pole, ahsahka$direction), 
                  colors=c(hues(ahsahka$ne), hues(ahsahka$ne)),
                  shapes=c(replicate(length(ahsahka$ne), "c"), replicate(length(ahsahka$ne), "s")),
                  curves=list(ahsahkaFols, ahsahkaLins))

# Both of those curves are small circles about a common axis. To get the axis 
# we must do a bit of math (sorry).
ahsahkaVort <- rotVectorFromAntisymmetric(t(ahsahkaRegr$b) %*% ahsahkaRegr$m %*% ahsahkaRegr$b)
geoTrendPlungeDegFromCartesian(rayNormalized(ahsahkaVort))

# Here is the amount of rotation about that axis, in degrees, per m of 
# northeasting. Or you can multiply it by 1000 to get the amount of rotation 
# per km of northeasting.
sqrt(dot(ahsahkaVort, ahsahkaVort)) / degree
sqrt(dot(ahsahkaVort, ahsahkaVort)) / degree * 1000

# We can also use the regression to predict missing values. For example, we 
# don't have any data at easting 560000 m, northing 5160000 m.
plot(x=ahsahka$easting, y=ahsahka$northing, xlab="easting (m)", ylab="northing (m)", col=hues(ahsahka$ne))

# Here's the predicted strike-dip of foliation and trend-plunge of lineation.
ahsahkaPred <- ahsahkaRegr$prediction((560000 + 5160000) / sqrt(2))
geoStrikeDipDegFromCartesian(ahsahkaPred[1,])
geoTrendPlungeDegFromCartesian(ahsahkaPred[2,])



### ANOTHER IDEA, POSSIBLY BAD ###

# After careful inspection of the plots above, we suspect that we might get 
# better results if we transformed the northeasting in certain ways. The 
# following transformation effectively exaggerates the differences between 
# northeastings near the shear zone. It will let the geodesic regression 
# produce big changes in orientation near the shear zone. (There's no other 
# geological motivation for this tranformation. It's just math.)
f <- function(ne) {-sin(2 * pi * (ne - min(ahsahka$ne)) / (max(ahsahka$ne) - min(ahsahka$ne)))}
ahsahkaXs <- sapply(ahsahka$ne, f)
plot(x=ahsahka$ne, y=ahsahkaXs, xlab="northeasting (m)", ylab="transformed northeasting")

# When we color the orientations by transformed northeasting, the rainbow 
# pattern is stronger than it was in the earlier plots.
oriEqualVolumePlot(ahsahka$rotation, group=oriLineInPlaneGroup, colors=hues(ahsahkaXs))

# Here's the geodesic regression. Check that the error is 0 and the 
# minEigenvalue is positive. The R^2 is better than it was earlier.
ahsahkaRegr <- oriRescaledGeodesicRegression(ahsahkaXs, ahsahka$rotation, group=oriLineInPlaneGroup, numSteps=1000)
ahsahkaRegr$error
ahsahkaRegr$minEigenvalue
ahsahkaRegr$rSquared

# Here is the prediction curve.
ahsahkaCurve <- lapply(seq(from=min(ahsahkaXs), to=max(ahsahkaXs), length.out=100), ahsahkaRegr$prediction)
oriEqualVolumePlot(ahsahka$rotation, group=oriLineInPlaneGroup, colors=hues(ahsahkaXs), curves=list(ahsahkaCurve), curveWidth=3)

# Here's the same plot, but with additional curves connecting each datum to 
# where the curve predicts it should be.
ahsahkaPreds <- lapply(ahsahkaXs, ahsahkaRegr$prediction)
ahsahkaCurves <- thread(function(pred, obs) rotGeodesicPoints(pred, oriNearestRepresentative(obs, pred, group=oriLineInPlaneGroup)),
                        ahsahkaPreds, ahsahka$rotation)
oriEqualVolumePlot(ahsahka$rotation, group=oriLineInPlaneGroup, colors=hues(ahsahkaXs), 
                   curves=c(ahsahkaCurves, list(ahsahkaCurve)), curveWidth=3)

# Here's the same plot, but emphasizing the predictions instead of the data.
oriEqualVolumePlot(ahsahkaPreds, group=oriLineInPlaneGroup, colors=hues(ahsahkaXs), 
                   curves=c(ahsahkaCurves, list(ahsahkaCurve)), curveWidth=3)

# Here's a plot of foliation poles (circles) and lineation directions 
# (squares) with the regression prediction curves superimposed.
ahsahkaFols <- lapply(ahsahkaCurve, function(r) r[1,])
ahsahkaLins <- lapply(ahsahkaCurve, function(r) r[2,])
lineEqualAreaPlot(c(ahsahka$pole, ahsahka$direction), 
                  colors=c(hues(ahsahkaXs), hues(ahsahkaXs)),
                  shapes=c(replicate(length(ahsahkaXs), "c"), replicate(length(ahsahkaXs), "s")),
                  curves=list(ahsahkaFols, ahsahkaLins))

# Both of those curves are small circles about a common axis.
ahsahkaVort <- rotVectorFromAntisymmetric(t(ahsahkaRegr$b) %*% ahsahkaRegr$m %*% ahsahkaRegr$b)
geoTrendPlungeDegFromCartesian(rayNormalized(ahsahkaVort))

# Here is the amount of rotation about that axis, in degrees, per unit of 
# transformed northeasting. But it might not mean much, because transformed 
# northeasting doesn't mean much.
sqrt(dot(ahsahkaVort, ahsahkaVort)) / degree

# Here's the permutation test, if you're interested.
ahsahkaRSqs <- permutedRSquareds(numPerms=100, oriRescaledGeodesicRegression, ahsahkaXs, 
                                 ahsahka$rotation, group=oriLineInPlaneGroup, numSteps=1000)
length(ahsahkaRSqs)
sum(ahsahkaRSqs > ahsahkaRegr$rSquared) / length(ahsahkaRSqs)



### CONCLUSION ###

# In geology, it is common for data to exhibit spatial or temporal patterns. 
# Regression is one tool for quantifying such patterns. But you need to make 
# sure that you're not over-detecting patterns that aren't really there.


