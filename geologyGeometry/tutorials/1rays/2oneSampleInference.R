


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

# In statistics, 'inference' is a process of extrapolating from a data set to 
# the larger population that it represents. In particular, the mean of the 
# data is our best estimate of the population mean, but there is some 
# uncertainty in that estimate. Confidence regions and hypothesis tests are 
# two closely related ways of studying that uncertainty. This tutorial 
# develops confidence regions and hypothesis tests for rays.



### PROBLEM ###

# Here is a set of paleomagnetic directions from Cyprus, compiled from various 
# papers (Bonhommet et al., 1988; Allerton, 1989; Abelson et al., 2002; Granot 
# et al., 2006; Scott et al., 2013).
pmag <- geoDataFromFile("data/cyprusSouthInPmag.csv")
rayEqualAreaPlot(pmag$direction)

# Suppose, for the sake of argument, that a certain proposal about the 
# geologic development of Cyprus predicts that paleomagnetic directions should 
# have trend 306 degrees and plunge 20 degrees. Here's that direction in red.
pmagNull <- geoCartesianFromTrendPlungeDeg(c(306, 20))
rayEqualAreaPlotTwo(pmag$direction, list(pmagNull), colorB="red")

# Even if the proposal is correct, the data won't all point in exactly that 
# predicted direction, if only because of random noise. But can the data be 
# explained using only random noise about that predicted direction? In other 
# words, could that direction be the mean of the population, from which the 
# data were drawn?

# It is worth noting that all of our inference methods assume that the data 
# are 'independent and identically distributed' (IID). That is, they all arise 
# from the same underlying population, and they constitute a 'random' sample 
# from that population (as opposed to a systematically biased sample, for 
# example).



### FISHER DISTRIBUTION ###

# The methods in this section assume that the data arise from a Fisher 
# distribution. In general, they make different approximations and hence 
# perform best for differing data sets. For this data set, however, they all 
# produce similar results.

# This first method is the 'default'. For small sample sizes (around n = 30) 
# it's a bit aggressive (small confidence regions, over-rejection of null 
# hypotheses). But for n = 43 it should work pretty well.
length(pmag$direction)

# Here's the 95% confidence region about the sample mean. Because the null 
# hypothesis is outside the 95% confidence region, we reject it at the 95% 
# confidence level.
pmagInf <- rayFisherConfidence(pmag$direction)
rayEqualAreaRadiusPlot(list(pmagInf$muHat, pmagNull), c(pmagInf$angle, 0), colors=c("black", "red"))

# Here's the radius of that confidence region --- an angle that is sometimes 
# called 'alpha_95'. And here's the angular distance from the center to the 
# null hypothesis. It's greater than alpha_95, confirming that the null 
# hypothesis is outside the region.
pmagInf$angle / degree
rayDistance(pmagInf$muHat, pmagNull) / degree

# The second method assumes large sample size, but n = 43 is again large 
# enough to assure good performance.
pmagInf <- rayFisherLargeSampleConfidence(pmag$direction)
rayEqualAreaRadiusPlot(list(pmagInf$muHat, pmagNull), c(pmagInf$angle, 0), colors=c("black", "red"))
pmagInf$angle / degree
rayDistance(pmagInf$muHat, pmagNull) / degree

# The third method isn't very sensitive to sample size, but instead works well 
# when the data are sufficiently concentrated --- such as kappa > 3. So we 
# check that kappa is at least 3 (yes, kappaHat ~ 4.5) and proceed as before.
pmagInf <- rayFisherTauxe(pmag$direction)
pmagInf$kappaHat
rayEqualAreaRadiusPlot(list(pmagInf$muHat, pmagNull), c(pmagInf$angle, 0), colors=c("black", "red"))
pmagInf$angle / degree
rayDistance(pmagInf$muHat, pmagNull) / degree



### BOOTSTRAPPING ###

# Bootstrapping is a simulation technique that makes no assumption about the 
# distribution from which the data arose. The idea is simple: We resample the 
# data with replacement, to create a mutated version of the data set, in which 
# some of the data are repeated and others are omitted. The mean of this 
# mutated data set is close to the original mean, but a little off. We 
# repeatedly resample and compute the mean, to build up a large set of these 
# mutated means. Bootstrapping theory (e.g., Efron and Tibshirani, 1993) says 
# that they quantify the uncertainty in the sample mean as an estimate of the 
# population mean.

# Here are 10,000 bootstrapped means, with the null hypothesis again in red.
pmagBoot <- rayBootstrapInference(pmag$direction, numBoots=10000, numPoints=100)
rayEqualAreaPlotTwo(pmagBoot$us, list(pmagNull), colorB="red")

# Here is an ellipse that contains 95% of the bootstrapped means. We take this 
# ellipse to be the 95% confidence region for the mean.
rayEqualAreaPlot(list(pmagBoot$center, pmagNull), curves=list(pmagBoot$points), colors=c("black", "red"))

# Because the null hypothesis is outside the 95% confidence region, we reject 
# the null hypothesis with 95% confidence. In fact, our bootstrapping 
# machinery provides an R function that produces a precise p-value. The fact 
# that p < 0.05 confirms that the null hypothesis is outside the region.
pmagBoot$pvalue(pmagNull)



### HOW TO REPORT AN INFERENCE ###

# When publishing inference results in a paper, do not merely report that you 
# reject or fail-to-reject the null hypothesis. Also report:
# 
# A. the specific method that you are using
# B. any important details --- for example, when bootstrapping, the number of 
#    bootstrap samples
# C. your chosen confidence level (usually 95%)
# D. the p-value (for methods that supply it)
# E. some idea of the effect size, such as a picture showing how far away from 
#    the confidence region the null hypothesis is.
# 
# The goal is to report enough information that the reader can determine 
# whether she agrees with your argument and reproduce your calculations if 
# necessary.



### OTHER OPTIONS ###

# Instead of using methods specific to the Fisher distribution, we could use 
# methods specific to the more flexible Kent distribution. Some tools are 
# available in the R package Directional. They have not yet been integrated 
# into our geologyGeometry library.

# Another option is Bayesian Markov chain Monte Carlo simulation. However, 
# this software library does not provide an easy, automated tool for doing 
# that technique on rays.



### CONCLUSION ###

# Bootstrapping is a pretty easy, reliable way to get at the population mean. 
# Fisher-based methods are also popular.


