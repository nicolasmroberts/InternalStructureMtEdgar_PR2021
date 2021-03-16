


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

# In general statistics, we frequently wish to extrapolate from our data set, 
# to make inferences about the wider population that it represents. The most 
# popular examples are confidence intervals and hypothesis tests for the 
# population mean. In this tutorial we study how to construct confidence 
# regions and perform hypothesis tests for orientational data.



### PROBLEM ###

# The data consist of 23 foliation-lineation pairs from the western Idaho 
# shear zone (Giorgis and Tikoff, 2004). Here are the 23 stations in map view.
wiszData <- geoDataFromFile("data/wiszFollins.tsv")
plot(x=wiszData$easting, y=wiszData$northing)

# The foliation-lineations have already been corrected for the Miocene-Recent 
# extension, by being rotated by the smallest rotation that takes a plane with 
# strike-dip (174, 25) degrees to horizontal. There's an outlier lineation, 
# which shows up as an outlier orientation.
lineEqualAreaPlotTwo(wiszData$pole, wiszData$direction, shapeB="s")
oriEqualAnglePlot(wiszData$rotation, group=oriLineInPlaneGroup)

# Giorgis and Tikoff (2004) model the western Idaho shear zone deformation as 
# homogeneous monoclinic transpression along a vertical, NS-striking shear 
# plane. This kind of deformation is controlled by two parameters: gamma 
# measures the amount of simple shear, and log k measures the amount of 
# shortening. Here are the foliation-lineation pairs that such a deformation 
# predicts, plotted in red.
monoPrediction <- function(gamma, logK) {
  f <- defMonoclinicPGT(gamma, logK)
  f <- rotMatrixAboutZ(pi / 2) %*% f %*% rotMatrixAboutZ(-pi / 2)
  eig <- eigen(f %*% t(f), symmetric=TRUE)
  fol <- eig$vectors[,3]
  lin <- eig$vectors[,1]
  rotProjectedMatrix(rbind(fol, lin, cross(fol, lin)))
}
monoPredictions <- list()
for (gamma in seq(from=0, to=50, by=0.5))
  for (logK in seq(from=0, to=-10, by=-0.5))
    monoPredictions[[length(monoPredictions) + 1]] <- monoPrediction(gamma, logK)
monoFoliations <- lapply(monoPredictions, function(r) r[1,])
monoLineations <- lapply(monoPredictions, function(r) r[2,])
oriEqualAnglePlot(c(wiszData$rotation, monoPredictions), 
                  colors=c(replicate(length(wiszData$rotation), "white"), replicate(length(monoPredictions), "red")),
                  group=oriLineInPlaneGroup, simplePoints=TRUE)

# The observed foliation poles look compatible with (at least some of) the 
# predicted foliation poles.
lineEqualAreaPlotTwo(wiszData$pole, monoFoliations, colorB="red")

# The observed lineation directions look compatible with the shortening-
# dominated transpressions (vertical lineation) but not the simple-shear-
# dominated transpressions (horizontal lineation). So Giorgis and Tikoff 
# (2004) inferred that the transpression was shortening-dominated.
lineEqualAreaPlotTwo(wiszData$direction, monoLineations, colorB="red", shapeA="s", shapeB="s")

# However, the observed lineations are not exactly vertical, right? But then 
# real data are never going to match theory exactly. There is always some 
# noise or error in real data.

# But it seems to me that the lineations are not really even centered on 
# vertical. They seem to be clumped around another direction, plunging steeply 
# to the north, that might indicate a different deformation geometry. So we 
# are led naturally to a geologic question: Could the observed foliation-
# lineations have been produced by the proposed deformation? A statistical 
# interpretation of this question is: Could any of the predicted foliation-
# lineations be the mean of the population, from which the data were drawn?



### PREPARATION ###

# Before we attempt such an inferential problem, we need to check some things.

# First, are the data independent and identically distributed? Probably not, 
# really. Foliation-lineations that are geographically close are probably 
# correlated. So let's view them colored by various spatial coordinates, such 
# as northing. In this plot, there is not a clear rainbow, but neither are the 
# colors totally mixed up. So there is probably some geographic dependency.
oriEqualAnglePlot(wiszData$rotation, group=oriLineInPlaneGroup, colors=hues(wiszData$northing))

# Fortunately, the Giorgis and Tikoff (2004) model assumes that the 
# deformation is homogeneous and hence that the foliation-lineations are 
# constant over the region. Under that null hypothesis, any variability in the 
# data must be due to random noise alone. So, for testing that null 
# hypothesis, we may proceed under the assumption that the data are 
# independent and identically distributed.

# In theory, all inference methods work all the time. In practice, some of 
# them perform better than others for certain sample sizes and dispersions. To 
# choose our inference method well, we estimate the dispersion, using one of 
# the 'cheats' from the preceding tutorial. The matrix Fisher K-hat has 
# eigenvalues 33, 11, 0.000003.
wiszMean <- oriFrechetMean(wiszData$rotation, group=oriLineInPlaneGroup)
wiszMean <- oriNearestRepresentative(wiszMean, diag(c(1, 1, 1)), group=oriLineInPlaneGroup)
wiszData$rotation <- oriNearestRepresentatives(wiszData$rotation, wiszMean, group=oriLineInPlaneGroup)
wiszMLE <- rotFisherMLE(wiszData$rotation, numSteps=1000)
eigen(wiszMLE$kHat, symmetric=TRUE, only.values=TRUE)$values

# Based on the numerical experiments reported by Davis and Titus (2017), we 
# should be using Markov chain Monte Carlo simulation. That computation 
# requires compilation of the C part of our library, which some users may find 
# difficult. (It's done in one of our C tutorials.) So here we're going to use 
# a simple bootstrapping method, even though it tends to produce confidence 
# regions that are slightly too small for data sets like this one.



### CONFIDENCE REGION ###

# Bootstrapping has been explained in a few earlier tutorials, but in case 
# you skipped those, here it is again. We resample the data set with 
# replacement, to form a mutated data set, in which some of the data are 
# omitted and others repeated. Then we compute the mean of that data set. It's 
# like a mutated version of the original sample mean. We repeat this process 
# many times, to build up a large set of mutated means. The theory of 
# bootstrapping (e.g., Efron and Tibshirani, 1993) says that these mutated 
# means quantify the uncertainty in the original sample mean as an estimate of 
# the population mean.

# For the sake of time, we'll compute only 1,000 bootstrapped orientation 
# means, which takes about 10 seconds on my laptop. If you were seriously 
# doing this technique for a research problem, I would recommend 10,000 
# bootstraps.
wiszBoot <- oriBootstrapInference(wiszData$rotation, numBoots=1000, group=oriLineInPlaneGroup)
oriEqualAnglePlot(wiszBoot$bootstraps, group=oriLineInPlaneGroup, simplePoints=TRUE)

# Because the bootstrapped means form a single tight clump (shown four times) 
# in the equal-angle plot, we can safely ignore the symmetry and focus on one 
# copy of the results.
rotEqualAnglePlot(wiszBoot$bootstraps, simplePoints=TRUE)

# Here are those bootstrapped means in equal-area. The foliation poles are 
# nearly horizontal and the lineation directions are nearly vertical but NOT 
# vertical.
lineEqualAreaPlot(c(lapply(wiszBoot$bootstraps, function(r) r[1,]), lapply(wiszBoot$bootstraps, function(r) r[2,])), shapes=c("."))

# We construct an ellipsoidal region that contains 95% of the simulation 
# results. This is a 95% confidence region for the population mean.
rotEllipsoidPlot(wiszBoot$bootstraps, wiszBoot$center, wiszBoot$leftCovarInv, wiszBoot$q095^2, numNonAdapt=4, simplePoints=TRUE)



### HYPOTHESIS TESTS ###

# Remember the foliation-lineations predicted by the monoclinic 
# transpressions? This plot shows them missing the 95% confidence region.
monoPredictions <- oriNearestRepresentatives(monoPredictions, wiszBoot$center, group=oriLineInPlaneGroup)
rotEllipsoidPlot(monoPredictions, wiszBoot$center, wiszBoot$leftCovarInv, wiszBoot$q095^2, numNonAdapt=4, simplePoints=TRUE, colors="red")

# In other words, if we did a hypothesis test with any one of these 
# predictions as the hypothesized mean, then that hypothesis would be rejected 
# with a p-value less than 0.05.  More precisely, here is the range of 
# p-values attained. (The particular numbers you see will depend on exactly 
# how your bootstrapping went.)
range(sapply(monoPredictions, wiszBoot$pvalueOri))

# So we reject the entire proposed class of deformations as an explanation for 
# these data.



### THE LANDSCAPE OF TRANSPRESSIONS ###

# We can picture the space of transpressions as a landscape with gamma running 
# east and logK running north. This first plot shows the 'critical curve' 
# where the finite strain ellipsoid is an oblate spheroid. Transpressions 
# below the curve are shortening-dominated and have vertical lineation. 
# Transpressions above the curve are simple-shear-dominated and have 
# horizontal lineation.
wiszLogKs <- seq(from=-2.4, to=0, by=0.01)
plot(x=sapply(wiszLogKs, defMonoclinicCriticalGammaFromLogK), y=wiszLogKs,
     pch=".", xlab="gamma", ylab="log k", main="transpressions")

# Maybe you think that we're picking on Giorgis and Tikoff (2004) too much. 
# Well, here are some results from Davis and Giorgis (2014). In that paper, we 
# tried to infer deformation from the rotation of rigid spheroidal clasts. 
# Each dot in this figure is a best-fit transpression for a perturbed version 
# of the clast data. (You might call this approach 'fake bootstrapping'.)
wisz2014Results <- geoDataFromFile("data/wiszGammaLogK.tsv")
points(x=wisz2014Results$gamma, y=wisz2014Results$logK)

# Our model did not incorporate the foliation-lineation data, but we knew that 
# the transpressions above the critical curve should be excluded for mis-
# predicting lineation. What we now know is that the transpressions below the 
# critical curve should also be excluded. The entire model is inappropriate.



### WHAT'S THE RIGHT ANSWER? ###

# When a statistical hypothesis test rejects its null hypothesis, it does not 
# give you a reason. In the Giorgis and Tikoff (2004) and Davis and Giorgis 
# (2014) cases, the null hypothesis is a fairly complicated 'narrative' about 
# how rocks in the western Idaho shear zone have behaved:
# * They were deformed by homogeneous transpression.
# * The shear plane was vertical and NS-striking.
# * The transport direction was horizontal.
# * Foliation and lineation align with the axes of the finite strain ellipsoid.
# * The observed foliation-lineations were independent and identically 
#   distributed (which is related to the homogeneity of the deformation).
# * There were no systematic measurement errors (compass declination set 
#   incorrectly, etc.).
# * And probably more...
# The hypothesis test does not tell you which part of this narrative is wrong. 
# Perhaps the transpression was heterogeneous (Robin and Cruden, 1994). 
# Perhaps it was triclinic (Lin et al., 1998; Jones and Holdsworth, 1998). 
# Perhaps it was oriented differently in space. Perhaps the foliation-
# lineation do not relate to finite strain exactly as we've assumed --- 
# especially near the critical curve of monoclinic transpression.

# Statistics doesn't tell you what the correct model is. The geologist must 
# use her creativity and intuition to invent a model that makes geologic 
# sense. She should then test whether the model is consistent with data. 
# Statistics is merely the 'calculator' for carrying out this consistency 
# check. If the model fails the check, then it should be discarded or altered. 
# If the model passes, then it can be accepted as provisionally true, with the 
# understanding that it may still be falsified in the future.



### CONCLUSION ###

# Confidence intervals and hypothesis tests are tools to help us weed out 
# models that are not compatible with the data.


